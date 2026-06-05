use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::{
    ExecutionContext, IntersectArgs, IntersectMode, log_verbose,
    vcf::{
        IndexRecord, OutputIndex, VariantKey, VcfRecord, open_text_reader, parse_header_samples,
        variant_keys_from_fields,
    },
};

#[derive(Default)]
struct IntersectMetrics {
    left_records: usize,
    right_records: usize,
    output_records: usize,
}

pub(crate) fn run(args: IntersectArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    match resolve_intersect_input(&args)? {
        IntersectInput::TwoFile { left, right } => {
            let left_record_count = count_records(&left)?;
            let right_record_count = count_records(&right)?;
            let right_keys = load_variant_keys(&right)?;
            let left_keys = if matches!(args.mode, IntersectMode::RightOnly) {
                load_variant_keys(&left)?
            } else {
                HashSet::new()
            };
            let mut metrics = match args.mode {
                IntersectMode::Shared => write_selected_records(
                    &left,
                    &args.output,
                    args.index_type,
                    &right_keys,
                    MatchSense::Present,
                    "both",
                )?,
                IntersectMode::LeftOnly => write_selected_records(
                    &left,
                    &args.output,
                    args.index_type,
                    &right_keys,
                    MatchSense::Absent,
                    "left",
                )?,
                IntersectMode::RightOnly => write_selected_records(
                    &right,
                    &args.output,
                    args.index_type,
                    &left_keys,
                    MatchSense::Absent,
                    "right",
                )?,
                IntersectMode::AllShared | IntersectMode::AnyShared | IntersectMode::SetDiff => {
                    unreachable!("multi-set modes are rejected before two-file execution")
                }
            };
            metrics.left_records = left_record_count;
            metrics.right_records = right_record_count;

            log_verbose(
                ctx,
                format!(
                    "intersect stage=done mode={:?} left_records={} right_records={} output_records={} output={} elapsed={:.2?}",
                    args.mode,
                    metrics.left_records,
                    metrics.right_records,
                    metrics.output_records,
                    args.output.display(),
                    started.elapsed()
                ),
            );
        }
        IntersectInput::OneFile {
            input,
            left_samples,
            right_samples,
        } => {
            let metrics = write_sample_group_records(
                &input,
                &args.output,
                args.index_type,
                args.mode,
                &left_samples,
                &right_samples,
            )?;
            log_verbose(
                ctx,
                format!(
                    "intersect stage=done mode={:?} input_records={} output_records={} output={} elapsed={:.2?}",
                    args.mode,
                    metrics.left_records,
                    metrics.output_records,
                    args.output.display(),
                    started.elapsed()
                ),
            );
        }
        IntersectInput::MultiSet {
            sets,
            emit_set,
            set_diff,
        } => {
            let metrics = write_multi_set_records(
                &sets,
                &emit_set,
                &args.output,
                args.index_type,
                args.mode,
                set_diff.as_ref(),
            )?;
            log_verbose(
                ctx,
                format!(
                    "intersect stage=done mode={:?} sets={} input_records={} output_records={} output={} elapsed={:.2?}",
                    args.mode,
                    sets.len(),
                    metrics.left_records,
                    metrics.output_records,
                    args.output.display(),
                    started.elapsed()
                ),
            );
        }
    }
    Ok(())
}

enum IntersectInput {
    TwoFile {
        left: PathBuf,
        right: PathBuf,
    },
    OneFile {
        input: PathBuf,
        left_samples: Vec<String>,
        right_samples: Vec<String>,
    },
    MultiSet {
        sets: Vec<NamedSetSpec>,
        emit_set: String,
        set_diff: Option<SetDiffExpr>,
    },
}

#[derive(Clone)]
struct NamedSetSpec {
    name: String,
    path: PathBuf,
    samples: Vec<String>,
}

#[derive(Clone)]
struct SetDiffExpr {
    left: String,
    right: String,
}

fn resolve_intersect_input(args: &IntersectArgs) -> Result<IntersectInput> {
    if !args.sets.is_empty() || args.set_manifest.is_some() {
        if args.input.is_some() || args.left.is_some() || args.right.is_some() {
            bail!("use either multi-set inputs or --input/--left/--right, not both");
        }
        if args.left_samples.is_some() || args.right_samples.is_some() {
            bail!("--left-samples and --right-samples require --input one-file mode");
        }
        let sets = collect_named_sets(&args.sets, args.set_manifest.as_deref())?;
        if sets.len() < 2 {
            bail!("multi-set mode requires at least two sets");
        }
        validate_unique_set_names(&sets)?;
        let set_diff = if matches!(args.mode, IntersectMode::SetDiff) {
            let expr = args
                .set_expr
                .as_deref()
                .context("set-diff mode requires a SET_EXPR such as A-B")?;
            Some(parse_set_diff_expr(expr)?)
        } else {
            if args.set_expr.is_some() {
                bail!("SET_EXPR is only supported with --mode set-diff");
            }
            None
        };
        let emit_set = args
            .emit_set
            .clone()
            .or_else(|| set_diff.as_ref().map(|expr| expr.left.clone()))
            .unwrap_or_else(|| sets[0].name.clone());
        if !sets.iter().any(|set| set.name == emit_set) {
            bail!("--emit-set names an unknown set: {emit_set}");
        }
        if let Some(expr) = &set_diff {
            validate_set_name_exists(&sets, &expr.left)?;
            validate_set_name_exists(&sets, &expr.right)?;
        }
        if matches!(
            args.mode,
            IntersectMode::Shared | IntersectMode::LeftOnly | IntersectMode::RightOnly
        ) {
            bail!("multi-set mode requires --mode all-shared, any-shared, or set-diff");
        }
        return Ok(IntersectInput::MultiSet {
            sets,
            emit_set,
            set_diff,
        });
    }

    if args.emit_set.is_some() || args.set_expr.is_some() {
        bail!("--emit-set and SET_EXPR require multi-set mode");
    }

    match (&args.input, &args.left, &args.right) {
        (Some(input), None, None) => {
            if matches!(
                args.mode,
                IntersectMode::AllShared | IntersectMode::AnyShared | IntersectMode::SetDiff
            ) {
                bail!("one-file mode requires --mode shared, left-only, or right-only");
            }
            let left_samples = args
                .left_samples
                .as_deref()
                .context("--left-samples is required with --input")
                .and_then(|value| parse_sample_list(value, "--left-samples"))?;
            let right_samples = args
                .right_samples
                .as_deref()
                .context("--right-samples is required with --input")
                .and_then(|value| parse_sample_list(value, "--right-samples"))?;
            Ok(IntersectInput::OneFile {
                input: input.clone(),
                left_samples,
                right_samples,
            })
        }
        (None, Some(left), Some(right)) => {
            if matches!(
                args.mode,
                IntersectMode::AllShared | IntersectMode::AnyShared | IntersectMode::SetDiff
            ) {
                bail!("two-file mode requires --mode shared, left-only, or right-only");
            }
            if args.left_samples.is_some() || args.right_samples.is_some() {
                bail!("--left-samples and --right-samples require --input one-file mode");
            }
            Ok(IntersectInput::TwoFile {
                left: left.clone(),
                right: right.clone(),
            })
        }
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
            bail!("use either --input or --left/--right, not both")
        }
        (None, _, _) => bail!("provide either --input or both --left and --right"),
    }
}

fn collect_named_sets(specs: &[String], manifest: Option<&Path>) -> Result<Vec<NamedSetSpec>> {
    let mut out = Vec::new();
    for spec in specs {
        out.push(parse_named_set_spec(spec)?);
    }
    if let Some(manifest) = manifest {
        let mut reader = open_text_reader(manifest)
            .with_context(|| format!("failed to open set manifest {}", manifest.display()))?;
        let mut line = String::new();
        let mut line_number = 0usize;
        while reader.read_line(&mut line)? != 0 {
            line_number += 1;
            let trimmed = line.trim_end_matches(['\r', '\n']);
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line.clear();
                continue;
            }
            out.push(parse_named_set_manifest_line(
                trimmed,
                manifest,
                line_number,
            )?);
            line.clear();
        }
    }
    Ok(out)
}

fn parse_named_set_spec(spec: &str) -> Result<NamedSetSpec> {
    let (name, rest) = spec
        .split_once('=')
        .with_context(|| format!("invalid --set {spec:?}; expected NAME=VCF[:SAMPLES]"))?;
    let (path, samples) = parse_set_path_and_samples(rest)?;
    build_named_set(name, path, samples)
}

fn parse_named_set_manifest_line(
    line: &str,
    manifest: &Path,
    line_number: usize,
) -> Result<NamedSetSpec> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if !(2..=3).contains(&fields.len()) {
        bail!(
            "invalid set manifest line {} in {}; expected NAME<TAB>VCF[<TAB>SAMPLES]",
            line_number,
            manifest.display()
        );
    }
    let samples = if fields.len() == 3 {
        parse_sample_list(fields[2], "manifest samples")?
    } else {
        Vec::new()
    };
    build_named_set(fields[0], PathBuf::from(fields[1]), samples)
}

fn parse_set_path_and_samples(rest: &str) -> Result<(PathBuf, Vec<String>)> {
    if let Some((path, samples)) = rest.rsplit_once(':') {
        return Ok((
            PathBuf::from(path),
            parse_sample_list(samples, "--set samples")?,
        ));
    }
    Ok((PathBuf::from(rest), Vec::new()))
}

fn build_named_set(name: &str, path: PathBuf, samples: Vec<String>) -> Result<NamedSetSpec> {
    validate_set_name(name)?;
    if path.as_os_str().is_empty() {
        bail!("set {name} has an empty VCF path");
    }
    Ok(NamedSetSpec {
        name: name.to_string(),
        path,
        samples,
    })
}

fn validate_set_name(name: &str) -> Result<()> {
    if name.is_empty() {
        bail!("set name cannot be empty");
    }
    if !name
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
    {
        bail!("set name {name:?} must contain only ASCII letters, numbers, '_', '-', or '.'");
    }
    Ok(())
}

fn validate_unique_set_names(sets: &[NamedSetSpec]) -> Result<()> {
    let mut names = HashSet::new();
    for set in sets {
        if !names.insert(set.name.as_str()) {
            bail!("duplicate set name: {}", set.name);
        }
    }
    Ok(())
}

fn validate_set_name_exists(sets: &[NamedSetSpec], name: &str) -> Result<()> {
    if !sets.iter().any(|set| set.name == name) {
        bail!("set expression names an unknown set: {name}");
    }
    Ok(())
}

fn parse_set_diff_expr(expr: &str) -> Result<SetDiffExpr> {
    let (left, right) = expr
        .split_once('-')
        .with_context(|| format!("invalid set-diff expression {expr:?}; expected A-B"))?;
    validate_set_name(left)?;
    validate_set_name(right)?;
    if left == right {
        bail!("set-diff expression must name two different sets");
    }
    Ok(SetDiffExpr {
        left: left.to_string(),
        right: right.to_string(),
    })
}

fn parse_sample_list(value: &str, flag: &str) -> Result<Vec<String>> {
    let samples = value
        .split(',')
        .map(str::trim)
        .filter(|sample| !sample.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    if samples.is_empty() {
        bail!("{flag} must name at least one sample");
    }
    Ok(samples)
}

#[derive(Clone, Copy)]
enum MatchSense {
    Present,
    Absent,
}

fn load_variant_keys(path: &Path) -> Result<HashSet<VariantKey>> {
    let mut reader =
        open_text_reader(path).with_context(|| format!("failed to open VCF {}", path.display()))?;
    let mut keys = HashSet::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            keys.extend(variant_keys_from_fields(&fields));
        }
        line.clear();
    }
    Ok(keys)
}

fn count_records(path: &Path) -> Result<usize> {
    let mut reader =
        open_text_reader(path).with_context(|| format!("failed to open VCF {}", path.display()))?;
    let mut count = 0usize;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            count += 1;
        }
        line.clear();
    }
    Ok(count)
}

fn write_selected_records(
    input: &Path,
    output: &Path,
    index_type: crate::IndexType,
    comparison_keys: &HashSet<VariantKey>,
    sense: MatchSense,
    membership: &str,
) -> Result<IntersectMetrics> {
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let mut reader = open_text_reader(input)
        .with_context(|| format!("failed to open input VCF {}", input.display()))?;
    let output_file = File::create(output)
        .with_context(|| format!("failed to create output {}", output.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    let mut output_index = OutputIndex::new(index_type);
    let mut metrics = IntersectMetrics::default();
    let mut wrote_header = false;
    let mut wrote_varlock_header = false;

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            if !wrote_varlock_header {
                write_varlock_intersect_headers(&mut writer)?;
                wrote_varlock_header = true;
            }
            writeln!(writer, "{trimmed}")?;
            wrote_header = true;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.left_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let keys = variant_keys_from_fields(&fields);
            let is_match = keys.iter().any(|key| comparison_keys.contains(key));
            let keep = match sense {
                MatchSense::Present => is_match,
                MatchSense::Absent => !is_match,
            };
            if keep {
                let annotated = add_membership_info(&fields, membership, None, None);
                let index_record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{annotated}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&index_record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }

    if !wrote_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

fn write_sample_group_records(
    input: &Path,
    output: &Path,
    index_type: crate::IndexType,
    mode: IntersectMode,
    left_samples: &[String],
    right_samples: &[String],
) -> Result<IntersectMetrics> {
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let mut reader = open_text_reader(input)
        .with_context(|| format!("failed to open input VCF {}", input.display()))?;
    let output_file = File::create(output)
        .with_context(|| format!("failed to create output {}", output.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    let mut output_index = OutputIndex::new(index_type);
    let mut metrics = IntersectMetrics::default();
    let mut wrote_header = false;
    let mut wrote_varlock_header = false;
    let mut sample_index = HashMap::new();

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            let (_, parsed_sample_index) = parse_header_samples(trimmed);
            validate_sample_group(left_samples, &parsed_sample_index, "--left-samples")?;
            validate_sample_group(right_samples, &parsed_sample_index, "--right-samples")?;
            sample_index = parsed_sample_index;
            if !wrote_varlock_header {
                write_varlock_intersect_headers(&mut writer)?;
                wrote_varlock_header = true;
            }
            writeln!(writer, "{trimmed}")?;
            wrote_header = true;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.left_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let record = VcfRecord::new(&fields, &sample_index);
            let left_support = group_alt_support_count(&record, left_samples)?;
            let right_support = group_alt_support_count(&record, right_samples)?;
            let membership = match mode {
                IntersectMode::Shared if left_support > 0 && right_support > 0 => Some("both"),
                IntersectMode::LeftOnly if left_support > 0 && right_support == 0 => Some("left"),
                IntersectMode::RightOnly if right_support > 0 && left_support == 0 => Some("right"),
                _ => None,
            };
            if let Some(membership) = membership {
                let annotated = add_membership_info(
                    &fields,
                    membership,
                    Some(left_support),
                    Some(right_support),
                );
                let index_record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{annotated}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&index_record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }

    if !wrote_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

struct NamedSetData {
    spec: NamedSetSpec,
    keys: HashSet<VariantKey>,
}

fn write_multi_set_records(
    sets: &[NamedSetSpec],
    emit_set: &str,
    output: &Path,
    index_type: crate::IndexType,
    mode: IntersectMode,
    set_diff: Option<&SetDiffExpr>,
) -> Result<IntersectMetrics> {
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let set_data = sets
        .iter()
        .map(|set| {
            load_named_set_keys(set).map(|keys| NamedSetData {
                spec: set.clone(),
                keys,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let emit = set_data
        .iter()
        .find(|set| set.spec.name == emit_set)
        .with_context(|| format!("unknown emit set: {emit_set}"))?;

    let mut reader = open_text_reader(&emit.spec.path)
        .with_context(|| format!("failed to open emit VCF {}", emit.spec.path.display()))?;
    let output_file = File::create(output)
        .with_context(|| format!("failed to create output {}", output.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    let mut output_index = OutputIndex::new(index_type);
    let mut metrics = IntersectMetrics::default();
    let mut wrote_header = false;
    let mut wrote_varlock_header = false;
    let mut sample_index = HashMap::new();

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            let (_, parsed_sample_index) = parse_header_samples(trimmed);
            if !emit.spec.samples.is_empty() {
                validate_sample_group(&emit.spec.samples, &parsed_sample_index, "--emit-set")?;
            }
            sample_index = parsed_sample_index;
            if !wrote_varlock_header {
                write_varlock_intersect_headers(&mut writer)?;
                wrote_varlock_header = true;
            }
            writeln!(writer, "{trimmed}")?;
            wrote_header = true;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.left_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let record_keys = record_keys_for_set_fields(&fields, &sample_index, &emit.spec)?;
            let matched_key = record_keys
                .iter()
                .find(|key| multi_set_key_selected(key, &set_data, mode, set_diff))
                .cloned();
            if let Some(key) = matched_key {
                let present_sets = set_data
                    .iter()
                    .filter(|set| set.keys.contains(&key))
                    .map(|set| set.spec.name.as_str())
                    .collect::<Vec<_>>();
                let annotated = add_multi_set_info(&fields, &present_sets);
                let index_record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{annotated}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&index_record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }

    if !wrote_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

fn load_named_set_keys(set: &NamedSetSpec) -> Result<HashSet<VariantKey>> {
    let mut reader = open_text_reader(&set.path)
        .with_context(|| format!("failed to open set VCF {}", set.path.display()))?;
    let mut keys = HashSet::new();
    let mut sample_index = HashMap::new();
    let mut saw_header = false;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            let (_, parsed_sample_index) = parse_header_samples(trimmed);
            if !set.samples.is_empty() {
                validate_sample_group(&set.samples, &parsed_sample_index, "--set samples")?;
            }
            sample_index = parsed_sample_index;
            saw_header = true;
        } else if !trimmed.is_empty() && !trimmed.starts_with('#') {
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            keys.extend(record_keys_for_set_fields(&fields, &sample_index, set)?);
        }
        line.clear();
    }
    if !saw_header {
        bail!(
            "set VCF is missing #CHROM header line: {}",
            set.path.display()
        );
    }
    Ok(keys)
}

fn record_keys_for_set_fields(
    fields: &[&str],
    sample_index: &HashMap<String, usize>,
    set: &NamedSetSpec,
) -> Result<Vec<VariantKey>> {
    if set.samples.is_empty() {
        return Ok(variant_keys_from_fields(fields));
    }

    let record = VcfRecord::new(fields, sample_index);
    let mut alt_indices = HashSet::new();
    for sample in &set.samples {
        alt_indices.extend(sample_alt_allele_indices(&record, sample)?);
    }
    Ok(variant_keys_for_alt_indices(fields, &alt_indices))
}

fn sample_alt_allele_indices(record: &VcfRecord<'_>, sample: &str) -> Result<HashSet<usize>> {
    let Some(gt) = record.sample_field(sample, "GT") else {
        return Ok(HashSet::new());
    };
    let mut out = HashSet::new();
    if gt == "." || gt == "./." || gt == ".|." || gt.is_empty() {
        return Ok(out);
    }
    for allele in gt.split(['/', '|']) {
        if allele == "." || allele.is_empty() {
            continue;
        }
        let allele_index = allele
            .parse::<usize>()
            .with_context(|| format!("invalid FORMAT/GT allele {allele:?} for sample {sample}"))?;
        if allele_index > 0 {
            out.insert(allele_index);
        }
    }
    Ok(out)
}

fn variant_keys_for_alt_indices(fields: &[&str], alt_indices: &HashSet<usize>) -> Vec<VariantKey> {
    let Some(chrom) = fields.first() else {
        return Vec::new();
    };
    let Some(pos) = fields.get(1) else {
        return Vec::new();
    };
    let Some(ref_allele) = fields.get(3) else {
        return Vec::new();
    };
    let Some(alt_field) = fields.get(4) else {
        return Vec::new();
    };
    if alt_field.is_empty() || *alt_field == "." {
        return Vec::new();
    }

    alt_field
        .split(',')
        .enumerate()
        .filter(|(i, alt)| alt_indices.contains(&(i + 1)) && !alt.is_empty() && *alt != ".")
        .map(|(_, alt)| VariantKey {
            chrom: (*chrom).to_string(),
            pos: (*pos).to_string(),
            ref_allele: (*ref_allele).to_string(),
            alt_allele: alt.to_string(),
        })
        .collect()
}

fn multi_set_key_selected(
    key: &VariantKey,
    sets: &[NamedSetData],
    mode: IntersectMode,
    set_diff: Option<&SetDiffExpr>,
) -> bool {
    match mode {
        IntersectMode::AllShared => sets.iter().all(|set| set.keys.contains(key)),
        IntersectMode::AnyShared => sets.iter().filter(|set| set.keys.contains(key)).count() >= 2,
        IntersectMode::SetDiff => {
            let Some(expr) = set_diff else {
                return false;
            };
            let left_present = sets
                .iter()
                .find(|set| set.spec.name == expr.left)
                .map(|set| set.keys.contains(key))
                .unwrap_or(false);
            let right_present = sets
                .iter()
                .find(|set| set.spec.name == expr.right)
                .map(|set| set.keys.contains(key))
                .unwrap_or(false);
            left_present && !right_present
        }
        IntersectMode::Shared | IntersectMode::LeftOnly | IntersectMode::RightOnly => false,
    }
}

fn validate_sample_group(
    samples: &[String],
    sample_index: &HashMap<String, usize>,
    flag: &str,
) -> Result<()> {
    for sample in samples {
        if !sample_index.contains_key(sample) {
            bail!("{flag} contains sample not found in VCF header: {sample}");
        }
    }
    Ok(())
}

fn group_alt_support_count(record: &VcfRecord<'_>, samples: &[String]) -> Result<usize> {
    let mut count = 0usize;
    for sample in samples {
        if sample_has_alt_gt(record, sample)? {
            count += 1;
        }
    }
    Ok(count)
}

fn sample_has_alt_gt(record: &VcfRecord<'_>, sample: &str) -> Result<bool> {
    let Some(gt) = record.sample_field(sample, "GT") else {
        return Ok(false);
    };
    if gt == "." || gt == "./." || gt == ".|." || gt.is_empty() {
        return Ok(false);
    }
    for allele in gt.split(['/', '|']) {
        if allele == "." || allele.is_empty() {
            continue;
        }
        let allele_index = allele
            .parse::<usize>()
            .with_context(|| format!("invalid FORMAT/GT allele {allele:?} for sample {sample}"))?;
        if allele_index > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn write_varlock_intersect_headers<W: Write>(writer: &mut W) -> Result<()> {
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_SET,Number=1,Type=String,Description=\"varlock intersect membership: left, right, or both\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_LEFT_SUPPORT,Number=1,Type=Integer,Description=\"Number of left sample-group genotypes with a non-reference FORMAT/GT allele\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_RIGHT_SUPPORT,Number=1,Type=Integer,Description=\"Number of right sample-group genotypes with a non-reference FORMAT/GT allele\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_SET_COUNT,Number=1,Type=Integer,Description=\"Number of named varlock intersect sets containing the emitted variant key\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_SETS,Number=.,Type=String,Description=\"Named varlock intersect sets containing the emitted variant key\">"
    )?;
    Ok(())
}

fn add_membership_info(
    fields: &[&str],
    membership: &str,
    left_support: Option<usize>,
    right_support: Option<usize>,
) -> String {
    let mut out = fields.to_vec();
    let mut info = if out[7] == "." || out[7].is_empty() {
        String::new()
    } else {
        out[7].to_string()
    };
    if !info.is_empty() {
        info.push(';');
    }
    info.push_str("VARLOCK_SET=");
    info.push_str(membership);
    if let Some(left_support) = left_support {
        info.push_str(";VARLOCK_LEFT_SUPPORT=");
        info.push_str(&left_support.to_string());
    }
    if let Some(right_support) = right_support {
        info.push_str(";VARLOCK_RIGHT_SUPPORT=");
        info.push_str(&right_support.to_string());
    }
    out[7] = &info;
    out.join("\t")
}

fn add_multi_set_info(fields: &[&str], present_sets: &[&str]) -> String {
    let mut out = fields.to_vec();
    let mut info = if out[7] == "." || out[7].is_empty() {
        String::new()
    } else {
        out[7].to_string()
    };
    if !info.is_empty() {
        info.push(';');
    }
    info.push_str("VARLOCK_SET=multi;VARLOCK_SET_COUNT=");
    info.push_str(&present_sets.len().to_string());
    info.push_str(";VARLOCK_SETS=");
    info.push_str(&present_sets.join(","));
    out[7] = &info;
    out.join("\t")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn intersect_shared_writes_left_records_with_membership() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let output = dir.path().join("shared.vcf.gz");
        write_vcf(
            &left,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.2\n",
        )?;
        write_vcf(&right, "chr1\t10\t.\tA\tC\t.\tPASS\t.\n")?;

        let args = IntersectArgs {
            left: Some(left),
            right: Some(right),
            input: None,
            left_samples: None,
            right_samples: None,
            sets: Vec::new(),
            set_manifest: None,
            emit_set: None,
            mode: IntersectMode::Shared,
            output: output.clone(),
            index_type: crate::IndexType::Csi,
            set_expr: None,
        };
        run(
            args,
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains("##INFO=<ID=VARLOCK_SET"));
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tVARLOCK_SET=both"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(dir.path().join("shared.vcf.gz.csi").exists());
        Ok(())
    }

    #[test]
    fn intersect_left_only_and_right_only_choose_output_side() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let left_only = dir.path().join("left_only.vcf.gz");
        let right_only = dir.path().join("right_only.vcf.gz");
        write_vcf(
            &left,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t11\t.\tA\tG\t.\tPASS\t.\n",
        )?;
        write_vcf(
            &right,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t12\t.\tA\tT\t.\tPASS\t.\n",
        )?;

        run(
            IntersectArgs {
                left: Some(left.clone()),
                right: Some(right.clone()),
                input: None,
                left_samples: None,
                right_samples: None,
                sets: Vec::new(),
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::LeftOnly,
                output: left_only.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        run(
            IntersectArgs {
                left: Some(left),
                right: Some(right),
                input: None,
                left_samples: None,
                right_samples: None,
                sets: Vec::new(),
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::RightOnly,
                output: right_only.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let left_text = read_bgzip(&left_only)?;
        assert!(left_text.contains("chr1\t11\t.\tA\tG\t.\tPASS\tVARLOCK_SET=left"));
        assert!(!left_text.contains("chr1\t12\t.\tA\tT"));

        let right_text = read_bgzip(&right_only)?;
        assert!(right_text.contains("chr1\t12\t.\tA\tT\t.\tPASS\tVARLOCK_SET=right"));
        assert!(!right_text.contains("chr1\t11\t.\tA\tG"));
        Ok(())
    }

    #[test]
    fn intersect_matches_any_alt_in_multialt_record() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let output = dir.path().join("shared.vcf.gz");
        write_vcf(&left, "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\n")?;
        write_vcf(&right, "chr1\t10\t.\tA\tG\t.\tPASS\t.\n")?;

        run(
            IntersectArgs {
                left: Some(left),
                right: Some(right),
                input: None,
                left_samples: None,
                right_samples: None,
                sets: Vec::new(),
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::Shared,
                output: output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains("chr1\t10\t.\tA\tC,G\t.\tPASS\tVARLOCK_SET=both"));
        Ok(())
    }

    #[test]
    fn intersect_one_file_shared_uses_sample_groups() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("cohort.vcf");
        let output = dir.path().join("shared.vcf.gz");
        write_multi_sample_vcf(
            &input,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n\
             chr1\t11\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n\
             chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\n",
        )?;

        run(
            IntersectArgs {
                left: None,
                right: None,
                input: Some(input),
                left_samples: Some("a,b".to_string()),
                right_samples: Some("c,d".to_string()),
                sets: Vec::new(),
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::Shared,
                output: output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains("##INFO=<ID=VARLOCK_LEFT_SUPPORT"));
        assert!(text.contains(
            "chr1\t10\t.\tA\tC\t.\tPASS\tVARLOCK_SET=both;VARLOCK_LEFT_SUPPORT=1;VARLOCK_RIGHT_SUPPORT=1\tGT\t0/1\t0/0\t0/1\t0/0"
        ));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(!text.contains("chr1\t12\t.\tA\tT"));
        Ok(())
    }

    #[test]
    fn intersect_one_file_left_only_treats_missing_and_absent_gt_as_absence() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("cohort.vcf");
        let output = dir.path().join("left_only.vcf.gz");
        write_multi_sample_vcf(
            &input,
            "chr1\t20\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t./.\t.|.\n\
             chr1\t21\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n\
             chr1\t22\t.\tA\tT\t.\tPASS\t.\tDP\t12\t10\t0\t0\n",
        )?;

        run(
            IntersectArgs {
                left: None,
                right: None,
                input: Some(input),
                left_samples: Some("a,b".to_string()),
                right_samples: Some("c,d".to_string()),
                sets: Vec::new(),
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::LeftOnly,
                output: output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains(
            "chr1\t20\t.\tA\tC\t.\tPASS\tVARLOCK_SET=left;VARLOCK_LEFT_SUPPORT=1;VARLOCK_RIGHT_SUPPORT=0\tGT\t0/1\t0/0\t./.\t.|."
        ));
        assert!(!text.contains("chr1\t21\t.\tA\tG"));
        assert!(!text.contains("chr1\t22\t.\tA\tT"));
        Ok(())
    }

    #[test]
    fn intersect_multi_set_all_shared_and_any_shared_use_named_sample_sets() -> Result<()> {
        let dir = tempdir()?;
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        let c = dir.path().join("c.vcf");
        let all_output = dir.path().join("all_shared.vcf.gz");
        let any_output = dir.path().join("any_shared.vcf.gz");
        write_multi_sample_vcf(
            &a,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n\
             chr1\t11\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n\
             chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n",
        )?;
        write_multi_sample_vcf(
            &b,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\n\
             chr1\t11\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\n",
        )?;
        write_multi_sample_vcf(
            &c,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\n\
             chr1\t13\t.\tA\tAAT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\n",
        )?;

        run(
            IntersectArgs {
                left: None,
                right: None,
                input: None,
                left_samples: None,
                right_samples: None,
                sets: vec![
                    format!("A={}:a", a.display()),
                    format!("B={}:b", b.display()),
                    format!("C={}:c", c.display()),
                ],
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::AllShared,
                output: all_output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        run(
            IntersectArgs {
                left: None,
                right: None,
                input: None,
                left_samples: None,
                right_samples: None,
                sets: vec![
                    format!("A={}:a", a.display()),
                    format!("B={}:b", b.display()),
                    format!("C={}:c", c.display()),
                ],
                set_manifest: None,
                emit_set: None,
                mode: IntersectMode::AnyShared,
                output: any_output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: None,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let all_text = read_bgzip(&all_output)?;
        assert!(all_text.contains("##INFO=<ID=VARLOCK_SETS"));
        assert!(all_text.contains(
            "chr1\t10\t.\tA\tC\t.\tPASS\tVARLOCK_SET=multi;VARLOCK_SET_COUNT=3;VARLOCK_SETS=A,B,C\tGT\t0/1\t0/0\t0/0\t0/0"
        ));
        assert!(!all_text.contains("chr1\t11\t.\tA\tG"));
        assert!(!all_text.contains("chr1\t12\t.\tA\tT"));

        let any_text = read_bgzip(&any_output)?;
        assert!(any_text.contains("chr1\t10\t.\tA\tC"));
        assert!(any_text.contains(
            "chr1\t11\t.\tA\tG\t.\tPASS\tVARLOCK_SET=multi;VARLOCK_SET_COUNT=2;VARLOCK_SETS=A,B\tGT\t0/1\t0/0\t0/0\t0/0"
        ));
        assert!(!any_text.contains("chr1\t12\t.\tA\tT"));
        Ok(())
    }

    #[test]
    fn intersect_multi_set_manifest_set_diff_uses_expression() -> Result<()> {
        let dir = tempdir()?;
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        let c = dir.path().join("c.vcf");
        let manifest = dir.path().join("sets.tsv");
        let output = dir.path().join("a_minus_b.vcf.gz");
        write_multi_sample_vcf(
            &a,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n\
             chr1\t11\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n",
        )?;
        write_multi_sample_vcf(
            &b,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\n",
        )?;
        write_multi_sample_vcf(
            &c,
            "chr1\t11\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\n",
        )?;
        std::fs::write(
            &manifest,
            format!(
                "A\t{}\ta\nB\t{}\tb\nC\t{}\tc\n",
                a.display(),
                b.display(),
                c.display()
            ),
        )?;

        run(
            IntersectArgs {
                left: None,
                right: None,
                input: None,
                left_samples: None,
                right_samples: None,
                sets: Vec::new(),
                set_manifest: Some(manifest),
                emit_set: None,
                mode: IntersectMode::SetDiff,
                output: output.clone(),
                index_type: crate::IndexType::Csi,
                set_expr: Some("A-B".to_string()),
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains(
            "chr1\t11\t.\tA\tG\t.\tPASS\tVARLOCK_SET=multi;VARLOCK_SET_COUNT=2;VARLOCK_SETS=A,C\tGT\t0/1\t0/0\t0/0\t0/0"
        ));
        assert!(!text.contains("chr1\t10\t.\tA\tC"));
        Ok(())
    }

    fn write_vcf(path: &Path, records: &str) -> Result<()> {
        std::fs::write(
            path,
            format!(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"
            ),
        )?;
        Ok(())
    }

    fn write_multi_sample_vcf(path: &Path, records: &str) -> Result<()> {
        std::fs::write(
            path,
            format!(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\n{records}"
            ),
        )?;
        Ok(())
    }

    fn read_bgzip(path: &Path) -> Result<String> {
        let mut reader = bgzf::io::Reader::new(File::open(path)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        Ok(text)
    }
}
