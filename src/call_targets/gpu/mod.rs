pub(crate) mod aggregate;
pub(crate) mod kernel;
pub(crate) mod runtime;
pub(crate) mod scan;

use std::{
    collections::{BTreeMap, HashMap},
    mem::size_of,
    sync::Arc,
    time::Instant,
};

use anyhow::{Context, Result, bail};

use crate::{CallTargetsGpuArgs, ExecutionContext, GpuBackend, call_targets, log_verbose};

use self::{
    aggregate::{build_chunk_plan, chunk_for_site, create_chunk_states, flush_chunk},
    kernel::create_kernel,
    runtime::{
        GpuSelector, auto_tuning_for_tier, classify_adapter, effective_matrix_budget,
        try_initialize_gpu,
    },
    scan::{
        CoveredObservation, CoveredScanEvent, ScanEvent, ScanWorkerParams,
        spawn_covered_scan_workers, spawn_scan_workers,
    },
};
use super::{
    derive_output,
    observation::{self, build_target_frontier_index, build_target_site_map},
    output,
    pileup::merge_counts,
    prepare_call_targets,
    types::{SiteCounts, SiteKey},
};

pub(crate) fn run(args: CallTargetsGpuArgs, ctx: &ExecutionContext) -> Result<()> {
    let label = "call_targets_gpu";
    let run_started = Instant::now();

    let selector = GpuSelector {
        list: args.gpu_list,
        indices: args.gpu_indices.clone(),
        name: args.gpu_name.clone(),
    };

    let Some(runtime) = try_initialize_gpu(backend_mask(args.gpu_backend), &selector, ctx.verbose)?
    else {
        if args.gpu_list {
            return Ok(());
        }
        if args.require_gpu {
            bail!("no compatible GPU adapter found");
        }
        log_verbose(
            ctx,
            format!("{label} no compatible GPU adapter found; falling back to CPU call-targets"),
        );
        return call_targets::run(args.call, ctx);
    };

    let prepared = prepare_call_targets(&args.call, ctx, label)?;
    let sample_count = prepared.sample_names.len();
    let output = args
        .call
        .output
        .clone()
        .unwrap_or_else(|| derive_output(&prepared.inputs, args.call.bamlist.as_deref()));
    log_verbose(ctx, format!("{label} output: {}", output.display()));

    let tier = classify_adapter(&runtime.adapter_info);
    let auto = auto_tuning_for_tier(tier);

    // Apply user overrides on top of tier defaults, then clamp to hardware limits.
    let requested_budget = args
        .matrix_budget_mib
        .map(|mib| mib.saturating_mul(1024 * 1024))
        .unwrap_or(auto.stream_matrix_budget_bytes);
    let requested_obs = args.max_obs_upload.unwrap_or(auto.max_obs_upload);

    let matrix_budget = effective_matrix_budget(&runtime.limits, requested_budget);
    let max_obs_upload = runtime::effective_max_obs_upload(
        &runtime.limits,
        requested_obs,
        size_of::<observation::Observation>(),
        kernel::WORKGROUP_SIZE as usize,
    );
    let flush_threshold = args.obs_flush_threshold.unwrap_or(max_obs_upload);

    if ctx.verbose > 0 {
        let auto_budget_mib = auto.stream_matrix_budget_bytes / (1024 * 1024);
        let eff_budget_mib = matrix_budget / (1024 * 1024);
        let budget_note = if args.matrix_budget_mib.is_some() {
            " (override)"
        } else {
            " (auto)"
        };
        let obs_note = if args.max_obs_upload.is_some() {
            " (override)"
        } else {
            " (auto)"
        };
        let flush_note = if args.obs_flush_threshold.is_some() {
            " (override)"
        } else {
            ""
        };
        eprintln!(
            "[{label}] adapter=\"{}\" tier={tier:?} \
             auto: matrix_budget={auto_budget_mib}MiB max_obs_upload={} \
             effective: matrix_budget={eff_budget_mib}MiB{budget_note} \
             max_obs_upload={max_obs_upload}{obs_note} \
             flush_threshold={flush_threshold}{flush_note}",
            runtime.adapter_info.name, auto.max_obs_upload,
        );
    }

    let kernel = create_kernel(&runtime, max_obs_upload)?;

    let all_counts = if args.call.targets.is_some() {
        run_static_target_gpu_path(
            label,
            ctx,
            &args,
            &runtime,
            &kernel,
            &prepared,
            sample_count,
            matrix_budget,
        )?
    } else {
        run_covered_gpu_path(
            label,
            ctx,
            &args,
            &runtime,
            &kernel,
            &prepared,
            sample_count,
            matrix_budget,
            flush_threshold,
        )?
    };

    output::write_call_targets_output(
        &args.call,
        ctx,
        label,
        &prepared.prepared_reference,
        &output,
        &prepared.ref_names,
        &prepared.sample_names,
        all_counts,
    )?;

    log_verbose(
        ctx,
        format!("{label} stage=done elapsed={:.2?}", run_started.elapsed()),
    );
    Ok(())
}

fn run_static_target_gpu_path(
    label: &str,
    ctx: &ExecutionContext,
    args: &CallTargetsGpuArgs,
    runtime: &runtime::GpuRuntime,
    kernel: &kernel::GpuAggregateKernel,
    prepared: &super::types::PreparedCallTargets,
    sample_count: usize,
    matrix_budget: usize,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    let target_started = Instant::now();
    let target_sites = Arc::new(build_target_site_map(&prepared.targets)?);
    if target_sites.is_empty() {
        bail!("no GPU target sites found");
    }
    let frontier_index = Arc::new(build_target_frontier_index(
        &target_sites.site_keys,
        prepared.ref_names.len(),
    )?);
    log_verbose(
        ctx,
        format!(
            "{label} stage=build_target_site_map sites={} elapsed={:.2?}",
            target_sites.len(),
            target_started.elapsed()
        ),
    );

    let chunk_plan = build_chunk_plan(target_sites.len(), sample_count, matrix_budget)?;
    let mut chunks = create_chunk_states(&chunk_plan, &runtime, sample_count, target_sites.len())?;
    log_verbose(
        ctx,
        format!(
            "{label} chunks={} max_sites_per_chunk={} matrix_budget={}",
            chunk_plan.total_chunks, chunk_plan.max_sites_per_chunk, matrix_budget
        ),
    );

    let scan_started = Instant::now();
    let scan_params = ScanWorkerParams {
        rg_to_sm: Arc::new(prepared.rg_to_sm.clone()),
        sample_index: Arc::new(prepared.sample_index.clone()),
        min_mapq: args.call.min_mapq,
        min_baseq: args.call.min_baseq,
        max_depth: args.call.max_depth,
        verbose: ctx.verbose,
    };
    let (scan_rx, scan_handles) = spawn_scan_workers(
        &prepared.inputs,
        scan_params,
        Arc::clone(&target_sites),
        frontier_index,
    );
    let mut done_inputs = 0usize;
    let mut skipped_flags = 0usize;
    let mut skipped_rg = 0usize;
    while done_inputs < prepared.inputs.len() {
        match scan_rx.recv() {
            Ok(ScanEvent::Batch { observations, .. }) => {
                for observation in observations {
                    let chunk_idx =
                        chunk_for_site(observation.site_idx, chunk_plan.max_sites_per_chunk);
                    let Some(chunk) = chunks.get_mut(chunk_idx) else {
                        bail!("scan emitted observation for invalid chunk {}", chunk_idx);
                    };
                    chunk.pending_obs.push(observation);
                }
            }
            Ok(ScanEvent::Progress { .. }) => {}
            Ok(ScanEvent::Done {
                skipped_rg: input_skipped_rg,
                skipped_flags: input_skipped_flags,
                ..
            }) => {
                done_inputs += 1;
                skipped_rg += input_skipped_rg;
                skipped_flags += input_skipped_flags;
            }
            Err(_) => break,
        }
    }
    for handle in scan_handles {
        handle
            .join()
            .map_err(|_| anyhow::anyhow!("GPU scan worker panicked"))??;
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=scan elapsed={:.2?} skipped_flags={} skipped_rg={}",
            scan_started.elapsed(),
            skipped_flags,
            skipped_rg
        ),
    );

    aggregate_chunks(
        label,
        ctx,
        args,
        runtime,
        kernel,
        &mut chunks,
        &target_sites.site_keys,
        sample_count,
    )
}

fn run_covered_gpu_path(
    label: &str,
    ctx: &ExecutionContext,
    args: &CallTargetsGpuArgs,
    runtime: &runtime::GpuRuntime,
    kernel: &kernel::GpuAggregateKernel,
    prepared: &super::types::PreparedCallTargets,
    sample_count: usize,
    matrix_budget: usize,
    flush_threshold: usize,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    log_verbose(
        ctx,
        format!(
            "{label} no --targets: streaming covered-site aggregation flush_threshold={}",
            flush_threshold
        ),
    );

    let scan_started = Instant::now();
    let scan_params = ScanWorkerParams {
        rg_to_sm: Arc::new(prepared.rg_to_sm.clone()),
        sample_index: Arc::new(prepared.sample_index.clone()),
        min_mapq: args.call.min_mapq,
        min_baseq: args.call.min_baseq,
        max_depth: args.call.max_depth,
        verbose: ctx.verbose,
    };
    let (scan_rx, scan_handles) = spawn_covered_scan_workers(&prepared.inputs, scan_params);

    let mut all_counts: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
    let mut pending: Vec<CoveredObservation> = Vec::new();
    let mut done_inputs = 0usize;
    let mut skipped_flags = 0usize;
    let mut skipped_rg = 0usize;
    let mut total_observations = 0u64;
    let mut flush_count = 0usize;

    while done_inputs < prepared.inputs.len() {
        match scan_rx.recv() {
            Ok(CoveredScanEvent::Batch { observations, .. }) => {
                total_observations += observations.len() as u64;
                pending.extend(observations);
                if pending.len() >= flush_threshold {
                    let batch = flush_covered_batch(
                        &pending,
                        kernel,
                        runtime,
                        sample_count,
                        args.call.max_depth,
                        matrix_budget,
                    )?;
                    merge_counts(&mut all_counts, batch, args.call.max_depth)?;
                    pending.clear();
                    flush_count += 1;
                    if ctx.verbose > 0 {
                        eprintln!(
                            "[{label}] covered flush #{flush_count}: obs_total={total_observations}"
                        );
                    }
                }
            }
            Ok(CoveredScanEvent::Done {
                skipped_rg: r,
                skipped_flags: f,
                ..
            }) => {
                done_inputs += 1;
                skipped_rg += r;
                skipped_flags += f;
            }
            Err(_) => break,
        }
    }
    for handle in scan_handles {
        handle
            .join()
            .map_err(|_| anyhow::anyhow!("GPU scan worker panicked"))??;
    }

    // Final flush for remaining observations.
    if !pending.is_empty() {
        let batch = flush_covered_batch(
            &pending,
            kernel,
            runtime,
            sample_count,
            args.call.max_depth,
            matrix_budget,
        )?;
        merge_counts(&mut all_counts, batch, args.call.max_depth)?;
        flush_count += 1;
    }

    log_verbose(
        ctx,
        format!(
            "{label} stage=scan+aggregate elapsed={:.2?} observations={} flush_count={} sites={} skipped_flags={} skipped_rg={}",
            scan_started.elapsed(),
            total_observations,
            flush_count,
            all_counts.len(),
            skipped_flags,
            skipped_rg
        ),
    );

    Ok(all_counts)
}

/// Remap, chunk, and GPU-aggregate one batch of covered observations.
/// Returns only sites with nonzero counts (sparse).
fn flush_covered_batch(
    covered: &[CoveredObservation],
    kernel: &kernel::GpuAggregateKernel,
    runtime: &runtime::GpuRuntime,
    sample_count: usize,
    max_depth: u32,
    matrix_budget: usize,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    let (site_keys, observations) = remap_covered_observations(covered.to_vec())?;
    if site_keys.is_empty() {
        return Ok(BTreeMap::new());
    }

    let chunk_plan = build_chunk_plan(site_keys.len(), sample_count, matrix_budget)?;
    let mut chunks = create_chunk_states(&chunk_plan, runtime, sample_count, site_keys.len())?;
    for obs in observations {
        let chunk_idx = chunk_for_site(obs.site_idx, chunk_plan.max_sites_per_chunk);
        let Some(chunk) = chunks.get_mut(chunk_idx) else {
            bail!("covered observation mapped to invalid chunk {}", chunk_idx);
        };
        chunk.pending_obs.push(obs);
    }

    let mut counts = BTreeMap::new();
    for chunk in &mut chunks {
        counts.extend(flush_chunk(
            chunk,
            kernel,
            runtime,
            &site_keys,
            sample_count,
            max_depth,
        )?);
    }
    Ok(counts)
}

fn remap_covered_observations(
    covered_observations: Vec<CoveredObservation>,
) -> Result<(Vec<SiteKey>, Vec<observation::Observation>)> {
    let mut site_keys = covered_observations
        .iter()
        .map(|observation| observation.site_key)
        .collect::<Vec<_>>();
    site_keys.sort_unstable();
    site_keys.dedup();

    let site_index = site_keys
        .iter()
        .enumerate()
        .map(|(idx, key)| {
            let idx = u32::try_from(idx).context("covered site count exceeds u32 range")?;
            Ok((*key, idx))
        })
        .collect::<Result<HashMap<_, _>>>()?;
    let observations = covered_observations
        .into_iter()
        .map(|covered| {
            let site_idx = *site_index
                .get(&covered.site_key)
                .context("covered site missing from dynamic site index")?;
            Ok(observation::Observation {
                site_idx,
                sample_base: covered.sample_base,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    Ok((site_keys, observations))
}

fn aggregate_chunks(
    label: &str,
    ctx: &ExecutionContext,
    args: &CallTargetsGpuArgs,
    runtime: &runtime::GpuRuntime,
    kernel: &kernel::GpuAggregateKernel,
    chunks: &mut [aggregate::GpuChunkState],
    site_keys: &[SiteKey],
    sample_count: usize,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    let aggregate_started = Instant::now();
    let mut all_counts: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
    for chunk in chunks {
        let chunk_counts = flush_chunk(
            chunk,
            kernel,
            runtime,
            site_keys,
            sample_count,
            args.call.max_depth,
        )?;
        all_counts.extend(chunk_counts);
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=aggregate elapsed={:.2?}",
            aggregate_started.elapsed()
        ),
    );

    Ok(all_counts)
}

fn backend_mask(backend: GpuBackend) -> wgpu::Backends {
    match backend {
        GpuBackend::All => wgpu::Backends::all(),
        GpuBackend::Vulkan => wgpu::Backends::VULKAN,
        GpuBackend::Metal => wgpu::Backends::METAL,
        GpuBackend::Dx12 => wgpu::Backends::DX12,
        GpuBackend::Gl => wgpu::Backends::GL,
    }
}
