mod output;
mod pileup;
mod reference;
mod samples;
mod targets;
mod types;

use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
    time::Instant,
};

use anyhow::{Context, Result, bail};
use noodles_bam::io::Reader;
use noodles_sam as sam;

use crate::{CallTargetsArgs, ExecutionContext, log_verbose};
use pileup::{InputResult, PileupSettings, ScanParams, merge_counts, process_input_bam};
use samples::{collect_samples, read_rg_map};
use targets::load_targets;
use types::{PreparedCallTargets, SiteCounts, SiteKey};

pub fn run(args: CallTargetsArgs, ctx: &ExecutionContext) -> Result<()> {
    let label = "call_targets";
    let prepared = prepare_call_targets(&args, ctx, label)?;
    let sample_count = prepared.sample_names.len();

    let scan_started = Instant::now();
    let mut all_counts: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
    for path in &prepared.inputs {
        let scan = ScanParams {
            rg_to_sm: &prepared.rg_to_sm,
            sample_index_map: &prepared.sample_index,
            min_mapq: args.min_mapq,
            verbose: ctx.verbose,
            pileup: PileupSettings {
                targets: &prepared.targets,
                sample_count,
                min_baseq: args.min_baseq,
                max_depth: args.max_depth,
            },
        };
        let InputResult {
            counts,
            skipped_flags,
            skipped_rg,
        } = process_input_bam(path, &scan)?;
        if ctx.verbose > 1 {
            eprintln!(
                "[{label}] scanned {} skipped_flags={} skipped_rg={}",
                path.display(),
                skipped_flags,
                skipped_rg
            );
        }
        merge_counts(&mut all_counts, counts, args.max_depth)?;
    }
    log_verbose(
        ctx,
        format!("{label} stage=scan elapsed={:.2?}", scan_started.elapsed()),
    );

    output::write_call_targets_output(
        &args,
        ctx,
        label,
        &prepared.prepared_reference,
        &prepared.ref_names,
        &prepared.sample_names,
        all_counts,
    )
}

pub(crate) fn prepare_call_targets(
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
) -> Result<PreparedCallTargets> {
    let stage_started = Instant::now();
    let prepared_reference =
        crate::fasta_prep::prepare_reference(&args.reference, ctx.verbose, ctx.threads)?;
    if prepared_reference != args.reference {
        log_verbose(
            ctx,
            format!(
                "{label} prepared reference: {} -> {}",
                args.reference.display(),
                prepared_reference.display()
            ),
        );
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=check_prepare_reference elapsed={:.2?}",
            stage_started.elapsed()
        ),
    );

    let stage_started = Instant::now();
    let inputs = resolve_inputs(args)?;
    if inputs.is_empty() {
        bail!("no input BAMs found");
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=resolve_inputs inputs={} elapsed={:.2?}",
            inputs.len(),
            stage_started.elapsed()
        ),
    );

    let stage_started = Instant::now();
    let (_header, ref_names, ref_name_to_id) = prepare_headers(&inputs)?;
    log_verbose(
        ctx,
        format!(
            "{label} stage=prepare_headers refs={} elapsed={:.2?}",
            ref_names.len(),
            stage_started.elapsed()
        ),
    );

    let stage_started = Instant::now();
    let targets = load_targets(&args.targets, &ref_name_to_id)?;
    if targets.by_ref.is_empty() {
        bail!("no target intervals match the BAM reference sequences");
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=load_targets target_refs={} elapsed={:.2?}",
            targets.by_ref.len(),
            stage_started.elapsed()
        ),
    );

    let stage_started = Instant::now();
    let rg_map = args
        .rg_map
        .as_ref()
        .map(|path| read_rg_map(path.as_path()))
        .transpose()?;
    let (sample_names, rg_to_sm) = collect_samples(&inputs, rg_map.as_deref())?;
    if sample_names.is_empty() {
        bail!("no samples found from RG/SM mapping");
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=collect_samples samples={} elapsed={:.2?}",
            sample_names.len(),
            stage_started.elapsed()
        ),
    );
    let sample_index = sample_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), i))
        .collect::<HashMap<_, _>>();

    log_verbose(ctx, format!("{label} inputs: {} BAMs", inputs.len()));
    log_verbose(ctx, format!("{label} samples: {:?}", sample_names));
    log_verbose(ctx, format!("{label} targets: {}", args.targets.display()));
    log_verbose(
        ctx,
        format!("{label} reference: {}", args.reference.display()),
    );

    Ok(PreparedCallTargets {
        prepared_reference,
        inputs,
        ref_names,
        targets,
        sample_names,
        rg_to_sm,
        sample_index,
    })
}

fn resolve_inputs(args: &CallTargetsArgs) -> Result<Vec<PathBuf>> {
    if !args.inputs.is_empty() {
        let mut paths = Vec::new();
        for input in &args.inputs {
            if input.is_dir() {
                collect_bams_recursive(input, &mut paths)?;
            } else {
                paths.push(input.clone());
            }
        }
        paths.sort();
        paths.dedup();
        return Ok(paths);
    }

    let bamlist = args
        .bamlist
        .as_ref()
        .context("bamlist required when no inputs are provided")?;
    let file = File::open(bamlist)
        .with_context(|| format!("failed to open bamlist {}", bamlist.display()))?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();
    for line in reader.lines() {
        let line = line.context("failed to read bamlist line")?;
        let path = line.trim();
        if path.is_empty() || path.starts_with('#') {
            continue;
        }
        paths.push(PathBuf::from(path));
    }
    paths.sort();
    paths.dedup();
    Ok(paths)
}

fn collect_bams_recursive(dir: &std::path::Path, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("failed to read directory {}", dir.display()))?
    {
        let entry = entry.context("failed to read directory entry")?;
        let path = entry.path();
        if path.is_dir() {
            collect_bams_recursive(&path, out)?;
        } else if path
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("bam"))
            .unwrap_or(false)
        {
            out.push(path);
        }
    }
    Ok(())
}

fn prepare_headers(
    paths: &[PathBuf],
) -> Result<(sam::Header, Vec<String>, HashMap<String, usize>)> {
    let mut reference_sequences = None;
    let mut header = None;
    for path in paths {
        let file = File::open(path)
            .with_context(|| format!("failed to open input BAM {}", path.display()))?;
        let mut reader = Reader::new(file);
        let input_header = reader
            .read_header()
            .with_context(|| format!("failed to read header for {}", path.display()))?;

        if let Some(ref refs) = reference_sequences {
            if refs != input_header.reference_sequences() {
                bail!(
                    "reference sequences for {} do not match the first input",
                    path.display()
                );
            }
        } else {
            reference_sequences = Some(input_header.reference_sequences().clone());
            header = Some(input_header.clone());
        }
    }

    let header = header.context("no inputs provided")?;
    let mut ref_names = Vec::new();
    let mut ref_name_to_id = HashMap::new();
    for (i, (name, _)) in header.reference_sequences().iter().enumerate() {
        ref_names.push(name.to_string());
        ref_name_to_id.insert(name.to_string(), i);
    }
    Ok((header, ref_names, ref_name_to_id))
}
