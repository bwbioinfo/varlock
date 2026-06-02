#![allow(dead_code)]

use std::collections::BTreeMap;

use anyhow::{Context, Result, bail};
use futures::{channel::oneshot, executor::block_on};

use super::{
    kernel::{self, GpuAggregateKernel},
    runtime::GpuRuntime,
};
use crate::call_targets::{
    observation::Observation,
    pileup::merge_sample_counts_with_cap,
    types::{SiteCounts, SiteKey},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ChunkPlan {
    pub(crate) max_sites_per_chunk: usize,
    pub(crate) total_chunks: usize,
    pub(crate) per_site_bytes: usize,
}

pub(crate) struct GpuChunkState {
    pub(crate) pending_obs: Vec<Observation>,
    pub(crate) counts_buffer: wgpu::Buffer,
    site_start: u32,
    site_end: u32,
}

pub(crate) fn build_chunk_plan(
    site_count: usize,
    sample_count: usize,
    matrix_budget_bytes: usize,
) -> Result<ChunkPlan> {
    if sample_count == 0 {
        bail!("sample_count must be nonzero for GPU chunk planning");
    }

    let per_site_bytes = sample_count
        .checked_mul(4)
        .and_then(|value| value.checked_mul(size_of::<u32>()))
        .context("GPU count matrix per-site byte size overflow")?;
    let max_sites_per_chunk = (matrix_budget_bytes / per_site_bytes).max(1);
    let total_chunks = site_count.div_ceil(max_sites_per_chunk);

    Ok(ChunkPlan {
        max_sites_per_chunk,
        total_chunks,
        per_site_bytes,
    })
}

pub(crate) fn chunk_for_site(site_idx: u32, max_sites_per_chunk: usize) -> usize {
    site_idx as usize / max_sites_per_chunk
}

pub(crate) fn chunk_site_range(
    chunk_idx: usize,
    max_sites_per_chunk: usize,
    site_count: usize,
) -> (usize, usize) {
    let start = chunk_idx.saturating_mul(max_sites_per_chunk);
    let end = start.saturating_add(max_sites_per_chunk).min(site_count);
    (start, end)
}

pub(crate) fn create_chunk_states(
    chunk_plan: &ChunkPlan,
    runtime: &GpuRuntime,
    sample_count: usize,
    site_count: usize,
) -> Result<Vec<GpuChunkState>> {
    let mut states = Vec::with_capacity(chunk_plan.total_chunks);
    for chunk_idx in 0..chunk_plan.total_chunks {
        states.push(create_chunk_state(
            chunk_idx,
            chunk_plan,
            runtime,
            sample_count,
            site_count,
        )?);
    }

    Ok(states)
}

pub(crate) fn create_chunk_state(
    chunk_idx: usize,
    chunk_plan: &ChunkPlan,
    runtime: &GpuRuntime,
    sample_count: usize,
    site_count: usize,
) -> Result<GpuChunkState> {
    let (site_start, site_end) =
        chunk_site_range(chunk_idx, chunk_plan.max_sites_per_chunk, site_count);
    let chunk_site_count = site_end.saturating_sub(site_start);
    let buffer_size = kernel::counts_buffer_size(chunk_site_count, sample_count)?;
    let counts_buffer = runtime.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("varlock.call_targets_gpu.aggregate.counts"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC
            | wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    if buffer_size > 0 {
        let zeroes = vec![0u8; usize::try_from(buffer_size).context("buffer too large")?];
        runtime.queue.write_buffer(&counts_buffer, 0, &zeroes);
    }

    Ok(GpuChunkState {
        pending_obs: Vec::new(),
        counts_buffer,
        site_start: u32::try_from(site_start).context("chunk site start exceeds u32 range")?,
        site_end: u32::try_from(site_end).context("chunk site end exceeds u32 range")?,
    })
}

pub(crate) fn flush_chunk(
    state: &mut GpuChunkState,
    kernel: &GpuAggregateKernel,
    runtime: &GpuRuntime,
    site_keys: &[SiteKey],
    sample_count: usize,
    max_depth: u32,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    let sample_count_u32 =
        u32::try_from(sample_count).context("sample count exceeds GPU parameter range")?;
    if !state.pending_obs.is_empty() {
        super::kernel::dispatch(
            kernel,
            runtime,
            &state.counts_buffer,
            &state.pending_obs,
            state.site_start,
            sample_count_u32,
        )?;
        state.pending_obs.clear();
    }

    let site_start = state.site_start as usize;
    let site_end = state.site_end as usize;
    let site_count = site_end.saturating_sub(site_start);
    let buffer_size = kernel::counts_buffer_size(site_count, sample_count)?;
    let counts = read_counts_buffer(runtime, &state.counts_buffer, buffer_size)?;
    materialize_counts(
        &counts,
        &site_keys[site_start..site_end],
        sample_count,
        max_depth,
    )
}

fn read_counts_buffer(
    runtime: &GpuRuntime,
    counts_buffer: &wgpu::Buffer,
    buffer_size: u64,
) -> Result<Vec<u32>> {
    if buffer_size == 0 {
        return Ok(Vec::new());
    }

    let staging = runtime.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("varlock.call_targets_gpu.aggregate.readback"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    let mut encoder = runtime
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("varlock.call_targets_gpu.aggregate.readback.encoder"),
        });
    encoder.copy_buffer_to_buffer(counts_buffer, 0, &staging, 0, buffer_size);
    runtime.queue.submit(Some(encoder.finish()));

    let slice = staging.slice(..);
    let (tx, rx) = oneshot::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    runtime
        .device
        .poll(wgpu::PollType::wait_indefinitely())
        .context("failed while waiting for GPU count readback")?;
    block_on(rx)
        .context("GPU count readback callback was dropped")?
        .context("failed to map GPU count readback buffer")?;

    let view = slice.get_mapped_range();
    let counts = bytemuck::cast_slice::<u8, u32>(&view).to_vec();
    drop(view);
    staging.unmap();

    Ok(counts)
}

fn materialize_counts(
    counts: &[u32],
    site_keys: &[SiteKey],
    sample_count: usize,
    max_depth: u32,
) -> Result<BTreeMap<SiteKey, SiteCounts>> {
    let expected = site_keys
        .len()
        .checked_mul(sample_count)
        .and_then(|value| value.checked_mul(4))
        .context("GPU count materialization size overflow")?;
    if counts.len() != expected {
        bail!(
            "GPU count readback length mismatch: got {}, expected {}",
            counts.len(),
            expected
        );
    }

    let mut out = BTreeMap::new();
    for (site_idx, site_key) in site_keys.iter().enumerate() {
        let mut per_sample = Vec::with_capacity(sample_count);
        let mut any_nonzero = false;
        for sample_idx in 0..sample_count {
            let offset = ((site_idx * sample_count) + sample_idx) * 4;
            let raw = [
                counts[offset],
                counts[offset + 1],
                counts[offset + 2],
                counts[offset + 3],
            ];
            let mut capped = [0u32; 4];
            merge_sample_counts_with_cap(&mut capped, raw, max_depth);
            any_nonzero |= capped.iter().any(|&value| value > 0);
            per_sample.push(capped);
        }

        if any_nonzero {
            out.insert(*site_key, SiteCounts { per_sample });
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(ref_id: usize, position: u32) -> SiteKey {
        SiteKey {
            reference_sequence_id: ref_id,
            position,
        }
    }

    #[test]
    fn build_chunk_plan_uses_sample_count_and_budget() -> Result<()> {
        let plan = build_chunk_plan(25, 2, 96)?;

        assert_eq!(
            plan,
            ChunkPlan {
                max_sites_per_chunk: 3,
                total_chunks: 9,
                per_site_bytes: 32,
            }
        );
        Ok(())
    }

    #[test]
    fn build_chunk_plan_rejects_zero_samples() {
        let err = build_chunk_plan(10, 0, 1024).unwrap_err();
        assert!(err.to_string().contains("sample_count"));
    }

    #[test]
    fn chunk_helpers_map_sites_to_ranges() {
        assert_eq!(chunk_for_site(0, 10), 0);
        assert_eq!(chunk_for_site(9, 10), 0);
        assert_eq!(chunk_for_site(10, 10), 1);
        assert_eq!(chunk_site_range(0, 10, 25), (0, 10));
        assert_eq!(chunk_site_range(2, 10, 25), (20, 25));
        assert_eq!(chunk_site_range(3, 10, 25), (30, 25));
    }

    #[test]
    fn materialize_counts_returns_sparse_sites_and_caps_depth() -> Result<()> {
        let site_keys = vec![site(0, 10), site(0, 11)];
        let counts = vec![
            0, 5, 5, 0, // site 0 sample 0 -> capped to depth 6
            0, 0, 0, 0, // site 0 sample 1
            0, 0, 0, 0, // site 1 sample 0
            0, 0, 0, 0, // site 1 sample 1
        ];

        let materialized = materialize_counts(&counts, &site_keys, 2, 6)?;

        assert_eq!(materialized.len(), 1);
        assert_eq!(materialized[&site(0, 10)].per_sample[0], [0, 3, 3, 0]);
        assert_eq!(materialized[&site(0, 10)].per_sample[1], [0, 0, 0, 0]);
        Ok(())
    }

    #[test]
    fn materialize_counts_errors_on_length_mismatch() {
        let err = materialize_counts(&[0, 1, 2], &[site(0, 1)], 1, 100).unwrap_err();
        assert!(err.to_string().contains("length mismatch"));
    }
}
