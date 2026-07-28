#[cfg(feature = "wgpu")]
pub(crate) mod gpu;
mod observation;
mod output;
mod pileup;
pub(crate) mod reference;
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
use pileup::{
    InputResult, PileupSettings, ScanParams, cap_counts, merge_counts, process_input_bam,
};
use samples::{collect_sample_resolution, read_rg_map};
use targets::load_targets;
use types::{
    Interval, PairedCallingConfig, PairedSampleRoles, PreparedCallTargets, SiteCounts, SiteKey,
    TargetIndex,
};

fn derive_output(resolved_inputs: &[PathBuf], bamlist: Option<&std::path::Path>) -> PathBuf {
    let source = resolved_inputs.first().map(|p| p.as_path()).or(bamlist);
    let stem = source
        .and_then(|p| p.file_stem())
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    PathBuf::from(format!("{stem}.vcf.gz"))
}

#[cfg(feature = "wgpu")]
pub fn run(args: CallTargetsArgs, ctx: &ExecutionContext) -> Result<()> {
    if args.gpu.cpu {
        return run_cpu(args, ctx);
    }

    gpu::run(crate::CallTargetsGpuArgs { call: args }, ctx)
}

#[cfg(not(feature = "wgpu"))]
pub fn run(args: CallTargetsArgs, ctx: &ExecutionContext) -> Result<()> {
    run_cpu(args, ctx)
}

pub(crate) fn run_cpu(args: CallTargetsArgs, ctx: &ExecutionContext) -> Result<()> {
    let label = "call_targets";
    let prepared = prepare_call_targets(&args, ctx, label)?;
    let sample_count = prepared.sample_names.len();

    let output = args
        .output
        .clone()
        .unwrap_or_else(|| derive_output(&prepared.inputs, args.bamlist.as_deref()));
    log_verbose(ctx, format!("{label} output: {}", output.display()));

    let scan_started = Instant::now();
    let mut all_counts: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
    for (path, sample_resolver) in prepared.inputs.iter().zip(&prepared.input_sample_resolvers) {
        let scan = ScanParams {
            sample_resolver,
            min_mapq: args.min_mapq,
            verbose: ctx.verbose,
            pileup: PileupSettings {
                targets: &prepared.targets,
                sample_count,
                min_baseq: args.min_baseq,
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
        merge_counts(&mut all_counts, counts)?;
    }
    // Apply the depth policy once after every input has contributed raw counts.
    cap_counts(&mut all_counts, args.max_depth);
    log_verbose(
        ctx,
        format!("{label} stage=scan elapsed={:.2?}", scan_started.elapsed()),
    );

    output::write_call_targets_output(
        &args,
        ctx,
        label,
        &prepared.prepared_reference,
        &output,
        &prepared.ref_names,
        &prepared.sample_names,
        prepared.paired.as_ref(),
        all_counts,
    )
}

pub(crate) fn prepare_call_targets(
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
) -> Result<PreparedCallTargets> {
    let stage_started = Instant::now();
    let reference = args
        .reference
        .as_ref()
        .context("--reference is required for call-targets")?;
    let prepared_reference =
        crate::fasta_prep::prepare_reference(reference, ctx.verbose, ctx.threads)?;
    if prepared_reference != *reference {
        log_verbose(
            ctx,
            format!(
                "{label} prepared reference: {} -> {}",
                reference.display(),
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
    let targets = match &args.targets {
        Some(path) => {
            let t = load_targets(path, &ref_name_to_id)?;
            if t.by_ref.is_empty() {
                bail!("no target intervals match the BAM reference sequences");
            }
            t
        }
        None => TargetIndex {
            by_ref: ref_name_to_id
                .values()
                .map(|&id| {
                    (
                        id,
                        vec![Interval {
                            start: 0,
                            end: u64::MAX,
                        }],
                    )
                })
                .collect(),
        },
    };
    log_verbose(
        ctx,
        format!(
            "{label} stage=load_targets target_refs={} elapsed={:.2?}",
            targets.by_ref.len(),
            stage_started.elapsed()
        ),
    );

    let stage_started = Instant::now();
    if args.rg_map.is_some() && args.split_by != crate::SplitBy::Sm {
        bail!("--rg-map is only supported with --split-by sm");
    }
    let rg_map = args
        .rg_map
        .as_ref()
        .map(|path| read_rg_map(path.as_path()))
        .transpose()?;
    let sample_resolution = collect_sample_resolution(&inputs, args.split_by, rg_map.as_deref())?;
    log_verbose(
        ctx,
        format!(
            "{label} stage=collect_samples samples={} elapsed={:.2?}",
            sample_resolution.sample_names.len(),
            stage_started.elapsed()
        ),
    );
    let paired = prepare_paired_calling(args, &sample_resolution.sample_index)?;

    log_verbose(ctx, format!("{label} inputs: {} BAMs", inputs.len()));
    log_verbose(
        ctx,
        format!(
            "{label} samples (--split-by {:?}): {:?}",
            args.split_by, sample_resolution.sample_names
        ),
    );
    for resolver in &sample_resolution.input_resolvers {
        log_verbose(
            ctx,
            format!(
                "{label} sample source: {}",
                resolver.describe(&sample_resolution.sample_names)
            ),
        );
    }
    if let Some(paired) = &paired {
        log_verbose(
            ctx,
            format!(
                "{label} paired: tumor={} normal={}",
                paired.roles.tumor, paired.roles.normal
            ),
        );
    }
    match &args.targets {
        Some(path) => log_verbose(ctx, format!("{label} targets: {}", path.display())),
        None => log_verbose(ctx, format!("{label} targets: all covered positions")),
    }
    log_verbose(ctx, format!("{label} reference: {}", reference.display()));

    Ok(PreparedCallTargets {
        prepared_reference,
        inputs,
        ref_names,
        targets,
        sample_names: sample_resolution.sample_names,
        input_sample_resolvers: sample_resolution.input_resolvers,
        paired,
    })
}

fn prepare_paired_calling(
    args: &CallTargetsArgs,
    sample_index: &HashMap<String, usize>,
) -> Result<Option<PairedCallingConfig>> {
    let Some(pair) = args.pair.as_deref() else {
        return Ok(None);
    };
    let roles = parse_pair_roles(pair)?;
    if roles.tumor == roles.normal {
        bail!("paired calling requires distinct tumor and normal samples");
    }
    let tumor_index = sample_index
        .get(&roles.tumor)
        .copied()
        .with_context(|| format!("paired tumor sample {:?} was not found", roles.tumor))?;
    let normal_index = sample_index
        .get(&roles.normal)
        .copied()
        .with_context(|| format!("paired normal sample {:?} was not found", roles.normal))?;

    Ok(Some(PairedCallingConfig {
        roles,
        tumor_index,
        normal_index,
        tumor_min_alt_count: args.tumor_min_alt_count,
        tumor_min_alt_fraction: args.tumor_min_alt_fraction,
        normal_max_alt_count: args.normal_max_alt_count,
        normal_max_alt_fraction: args.normal_max_alt_fraction,
        normal_min_depth: args.normal_min_depth,
    }))
}

fn parse_pair_roles(value: &str) -> Result<PairedSampleRoles> {
    let mut tumor = None;
    let mut normal = None;
    for part in value.split(',') {
        let (key, sample) = part
            .split_once('=')
            .with_context(|| format!("invalid --pair {value:?}; expected tumor=S,normal=S"))?;
        if sample.is_empty() {
            bail!("invalid --pair {value:?}; sample name cannot be empty");
        }
        match key {
            "tumor" => tumor = Some(sample.to_string()),
            "normal" => normal = Some(sample.to_string()),
            _ => bail!("invalid --pair role {key:?}; expected tumor or normal"),
        }
    }
    let tumor = tumor.with_context(|| format!("invalid --pair {value:?}; missing tumor=S"))?;
    let normal = normal.with_context(|| format!("invalid --pair {value:?}; missing normal=S"))?;
    Ok(PairedSampleRoles { tumor, normal })
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IndexType;

    fn call_targets_args_with_pair(pair: &str) -> CallTargetsArgs {
        CallTargetsArgs {
            inputs: Vec::new(),
            bamlist: None,
            reference: None,
            targets: None,
            output: None,
            rg_map: None,
            split_by: crate::SplitBy::Sm,
            #[cfg(feature = "wgpu")]
            gpu: crate::GpuArgs::default(),
            index_type: IndexType::Csi,
            min_mapq: 20,
            min_baseq: 20,
            min_alt_count: 1,
            min_alt_fraction: 0.0,
            pair: Some(pair.to_string()),
            tumor_min_alt_count: 2,
            tumor_min_alt_fraction: 0.1,
            normal_max_alt_count: 0,
            normal_max_alt_fraction: 0.01,
            normal_min_depth: Some(10),
            max_depth: 100_000,
        }
    }

    #[test]
    fn parse_pair_roles_accepts_tumor_and_normal() -> Result<()> {
        let roles = parse_pair_roles("tumor=Tumor,normal=Normal")?;
        assert_eq!(roles.tumor, "Tumor");
        assert_eq!(roles.normal, "Normal");
        Ok(())
    }

    #[test]
    fn parse_pair_roles_rejects_missing_role() {
        let err = parse_pair_roles("tumor=Tumor").unwrap_err();
        assert!(err.to_string().contains("missing normal"));
    }

    #[test]
    fn prepare_paired_calling_validates_samples() {
        let args = call_targets_args_with_pair("tumor=Tumor,normal=Normal");
        let sample_index = [
            ("Tumor".to_string(), 0usize),
            ("Normal".to_string(), 1usize),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();

        let paired = prepare_paired_calling(&args, &sample_index)
            .unwrap()
            .unwrap();
        assert_eq!(paired.tumor_index, 0);
        assert_eq!(paired.normal_index, 1);
        assert_eq!(paired.tumor_min_alt_count, 2);
        assert_eq!(paired.tumor_min_alt_fraction, 0.1);
        assert_eq!(paired.normal_max_alt_count, 0);
        assert_eq!(paired.normal_max_alt_fraction, 0.01);
        assert_eq!(paired.normal_min_depth, Some(10));
    }

    #[test]
    fn prepare_paired_calling_rejects_unknown_and_duplicate_samples() {
        let sample_index = [("Tumor".to_string(), 0usize)]
            .into_iter()
            .collect::<HashMap<_, _>>();

        let args = call_targets_args_with_pair("tumor=Tumor,normal=Missing");
        let err = prepare_paired_calling(&args, &sample_index).unwrap_err();
        assert!(err.to_string().contains("paired normal sample"));

        let args = call_targets_args_with_pair("tumor=Tumor,normal=Tumor");
        let err = prepare_paired_calling(&args, &sample_index).unwrap_err();
        assert!(err.to_string().contains("distinct"));
    }
}
