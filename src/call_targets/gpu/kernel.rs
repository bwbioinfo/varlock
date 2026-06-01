#![allow(dead_code)]

use anyhow::{Context, Result, bail};

use super::runtime::GpuRuntime;
use crate::call_targets::observation::Observation;

const WORKGROUP_SIZE: u32 = 256;

pub(crate) const AGGREGATE_SHADER: &str = r#"
struct Observation {
    site_idx: u32,
    sample_base: u32,
};

struct Params {
    sample_count: u32,
    obs_len: u32,
    site_start: u32,
    _pad: u32,
};

@group(0) @binding(0)
var<storage, read> observations: array<Observation>;
@group(0) @binding(1)
var<storage, read_write> counts: array<atomic<u32>>;
@group(0) @binding(2)
var<storage, read> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.obs_len) {
        return;
    }
    let obs = observations[i];
    let sample_idx = obs.sample_base >> 2u;
    let base_idx = obs.sample_base & 3u;
    let local_site = obs.site_idx - params.site_start;
    let out_idx = ((local_site * params.sample_count + sample_idx) * 4u) + base_idx;
    atomicAdd(&counts[out_idx], 1u);
}
"#;

pub(crate) struct GpuAggregateKernel {
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    obs_buffer: wgpu::Buffer,
    params_buffer: wgpu::Buffer,
    pub(crate) max_obs_upload: usize,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    sample_count: u32,
    obs_len: u32,
    site_start: u32,
    _pad: u32,
}

pub(crate) fn create_kernel(
    runtime: &GpuRuntime,
    max_obs_upload: usize,
) -> Result<GpuAggregateKernel> {
    let max_obs_upload = max_obs_upload.max(1);
    let shader = runtime
        .device
        .create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("varlock.call_targets_gpu.aggregate"),
            source: wgpu::ShaderSource::Wgsl(AGGREGATE_SHADER.into()),
        });
    let bind_group_layout =
        runtime
            .device
            .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("varlock.call_targets_gpu.aggregate.bind_group_layout"),
                entries: &[
                    storage_buffer_layout_entry(0, true),
                    storage_buffer_layout_entry(1, false),
                    storage_buffer_layout_entry(2, true),
                ],
            });
    let pipeline_layout = runtime
        .device
        .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("varlock.call_targets_gpu.aggregate.pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            immediate_size: 0,
        });
    let pipeline = runtime
        .device
        .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("varlock.call_targets_gpu.aggregate.pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

    let obs_buffer = runtime.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("varlock.call_targets_gpu.aggregate.observations"),
        size: (max_obs_upload * size_of::<Observation>()) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    let params_buffer = runtime.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("varlock.call_targets_gpu.aggregate.params"),
        size: size_of::<Params>() as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });

    Ok(GpuAggregateKernel {
        pipeline,
        bind_group_layout,
        obs_buffer,
        params_buffer,
        max_obs_upload,
    })
}

pub(crate) fn dispatch(
    kernel: &GpuAggregateKernel,
    runtime: &GpuRuntime,
    counts_buffer: &wgpu::Buffer,
    observations: &[Observation],
    site_start: u32,
    sample_count: u32,
) -> Result<()> {
    if sample_count == 0 {
        bail!("sample_count must be nonzero for GPU aggregation");
    }

    for chunk in observations.chunks(kernel.max_obs_upload) {
        let obs_len = u32::try_from(chunk.len()).context("observation batch exceeds u32 range")?;
        if obs_len == 0 {
            continue;
        }
        let params = Params {
            sample_count,
            obs_len,
            site_start,
            _pad: 0,
        };

        runtime
            .queue
            .write_buffer(&kernel.obs_buffer, 0, bytemuck::cast_slice(chunk));
        runtime
            .queue
            .write_buffer(&kernel.params_buffer, 0, bytemuck::bytes_of(&params));

        let bind_group = runtime
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("varlock.call_targets_gpu.aggregate.bind_group"),
                layout: &kernel.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: kernel.obs_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: counts_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: kernel.params_buffer.as_entire_binding(),
                    },
                ],
            });
        let mut encoder = runtime
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("varlock.call_targets_gpu.aggregate.encoder"),
            });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("varlock.call_targets_gpu.aggregate.pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&kernel.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(obs_len.div_ceil(WORKGROUP_SIZE), 1, 1);
        }

        runtime.queue.submit(Some(encoder.finish()));
        runtime
            .device
            .poll(wgpu::PollType::wait_indefinitely())
            .context("failed while waiting for GPU aggregate dispatch")?;
    }

    Ok(())
}

pub(crate) fn counts_buffer_size(site_count: usize, sample_count: usize) -> Result<u64> {
    let cells = site_count
        .checked_mul(sample_count)
        .and_then(|value| value.checked_mul(4))
        .context("GPU count matrix cell count overflow")?;
    let bytes = cells
        .checked_mul(size_of::<u32>())
        .context("GPU count matrix byte size overflow")?;
    u64::try_from(bytes).context("GPU count matrix byte size exceeds u64 range")
}

fn storage_buffer_layout_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_buffer_size_uses_sites_samples_and_four_bases() -> Result<()> {
        assert_eq!(counts_buffer_size(10, 3)?, 10 * 3 * 4 * 4);
        Ok(())
    }

    #[test]
    fn counts_buffer_size_detects_overflow() {
        let err = counts_buffer_size(usize::MAX, 2).unwrap_err();
        assert!(err.to_string().contains("overflow"));
    }
}
