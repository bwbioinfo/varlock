#![allow(dead_code)]

use anyhow::{Context, Result, bail};
use futures::executor::block_on;

pub(crate) struct GpuRuntime {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    pub(crate) adapter_info: wgpu::AdapterInfo,
    pub(crate) limits: wgpu::Limits,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct GpuSelector {
    pub(crate) list: bool,
    pub(crate) index: Option<usize>,
    pub(crate) name: Option<String>,
}

#[derive(Clone, Debug)]
pub(crate) struct GpuAdapterSummary {
    pub(crate) index: usize,
    pub(crate) info: wgpu::AdapterInfo,
    pub(crate) limits: wgpu::Limits,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum GpuTier {
    Datacenter,
    HighEndDiscrete,
    Discrete,
    IntegratedOrLaptop,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct GpuAutoTuning {
    pub(crate) stream_matrix_budget_bytes: usize,
    pub(crate) max_obs_upload: usize,
    pub(crate) ready_batch_obs_limit: usize,
}

pub(crate) fn try_initialize_gpu(
    backend_mask: wgpu::Backends,
    selector: &GpuSelector,
    verbose: u8,
) -> Result<Option<GpuRuntime>> {
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
        backends: backend_mask,
        ..Default::default()
    });

    if selector.list {
        let summaries = enumerate_adapter_summaries(&instance, backend_mask);
        print_adapter_summaries(&summaries);
        return Ok(None);
    }

    let adapter = if selector.index.is_some() || selector.name.is_some() {
        let adapters = block_on(instance.enumerate_adapters(backend_mask));
        let summaries = adapter_summaries(&adapters);
        let selected_idx = select_adapter_index(&summaries, selector)?;
        adapters
            .into_iter()
            .nth(selected_idx)
            .with_context(|| format!("selected GPU adapter index {selected_idx} disappeared"))?
    } else {
        match block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
        })) {
            Ok(adapter) => adapter,
            Err(_) => return Ok(None),
        }
    };

    let adapter_info = adapter.get_info();
    if verbose > 0 {
        eprintln!(
            "[call_targets_gpu] selected adapter: {} ({:?})",
            adapter_info.name, adapter_info.device_type
        );
    }

    let (device, queue) = block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: Some("varlock-call-targets-gpu"),
        required_features: wgpu::Features::empty(),
        required_limits: wgpu::Limits::default(),
        ..Default::default()
    }))
    .context("failed to create wgpu device")?;
    let limits = device.limits();

    Ok(Some(GpuRuntime {
        device,
        queue,
        adapter_info,
        limits,
    }))
}

fn enumerate_adapter_summaries(
    instance: &wgpu::Instance,
    backend_mask: wgpu::Backends,
) -> Vec<GpuAdapterSummary> {
    let adapters = block_on(instance.enumerate_adapters(backend_mask));
    adapter_summaries(&adapters)
}

fn adapter_summaries(adapters: &[wgpu::Adapter]) -> Vec<GpuAdapterSummary> {
    adapters
        .iter()
        .enumerate()
        .map(|(index, adapter)| GpuAdapterSummary {
            index,
            info: adapter.get_info(),
            limits: adapter.limits(),
        })
        .collect()
}

fn print_adapter_summaries(summaries: &[GpuAdapterSummary]) {
    if summaries.is_empty() {
        println!("No compatible GPU adapters found");
        return;
    }

    for summary in summaries {
        let info = &summary.info;
        let limits = &summary.limits;
        println!(
            "[{}] {} ({:?}, backend={:?}, vendor=0x{:04x}, device=0x{:04x}, pci_bus_id={})",
            summary.index,
            info.name,
            info.device_type,
            info.backend,
            info.vendor,
            info.device,
            empty_as_dash(&info.device_pci_bus_id)
        );
        println!(
            "    driver={} {} max_storage_binding={} max_buffer_size={} max_workgroups_x={}",
            empty_as_dash(&info.driver),
            empty_as_dash(&info.driver_info),
            limits.max_storage_buffer_binding_size,
            limits.max_buffer_size,
            limits.max_compute_workgroups_per_dimension
        );
    }
}

fn empty_as_dash(value: &str) -> &str {
    if value.is_empty() { "-" } else { value }
}

fn select_adapter_index(summaries: &[GpuAdapterSummary], selector: &GpuSelector) -> Result<usize> {
    if summaries.is_empty() {
        bail!("no compatible GPU adapters found");
    }
    if selector.index.is_some() && selector.name.is_some() {
        bail!("--gpu-index and --gpu-name cannot be used together");
    }
    if let Some(index) = selector.index {
        if summaries.iter().any(|summary| summary.index == index) {
            return Ok(index);
        }
        bail!(
            "GPU adapter index {} is out of range; {} adapter(s) available",
            index,
            summaries.len()
        );
    }
    if let Some(name) = &selector.name {
        let needle = name.to_ascii_lowercase();
        if let Some(summary) = summaries
            .iter()
            .find(|summary| summary.info.name.to_ascii_lowercase().contains(&needle))
        {
            return Ok(summary.index);
        }
        bail!("no GPU adapter name contains {:?}", name);
    }

    bail!("no GPU adapter selector was provided")
}

pub(crate) fn classify_adapter(info: &wgpu::AdapterInfo) -> GpuTier {
    let name = info.name.to_ascii_lowercase();

    if matches!(
        info.device_type,
        wgpu::DeviceType::IntegratedGpu
            | wgpu::DeviceType::VirtualGpu
            | wgpu::DeviceType::Cpu
            | wgpu::DeviceType::Other
    ) {
        return GpuTier::IntegratedOrLaptop;
    }

    if [
        "h200", "h100", "gh200", "a100", "a30", "a40", "b100", "b200", "v100", "p100", "l40",
        "l40s",
    ]
    .iter()
    .any(|needle| name.contains(needle))
    {
        return GpuTier::Datacenter;
    }

    if [
        "rtx 4090",
        "rtx 4080",
        "rtx 3090",
        "rtx 6000",
        "rtx a6000",
        "rtx a5000",
    ]
    .iter()
    .any(|needle| name.contains(needle))
    {
        return GpuTier::HighEndDiscrete;
    }

    if matches!(info.device_type, wgpu::DeviceType::DiscreteGpu) {
        return GpuTier::Discrete;
    }

    GpuTier::IntegratedOrLaptop
}

pub(crate) fn auto_tuning_for_tier(tier: GpuTier) -> GpuAutoTuning {
    match tier {
        GpuTier::Datacenter => GpuAutoTuning {
            stream_matrix_budget_bytes: 1024 * 1024 * 1024,
            max_obs_upload: 64_000_000,
            ready_batch_obs_limit: 750_000,
        },
        GpuTier::HighEndDiscrete => GpuAutoTuning {
            stream_matrix_budget_bytes: 512 * 1024 * 1024,
            max_obs_upload: 32_000_000,
            ready_batch_obs_limit: 350_000,
        },
        GpuTier::Discrete => GpuAutoTuning {
            stream_matrix_budget_bytes: 256 * 1024 * 1024,
            max_obs_upload: 24_000_000,
            ready_batch_obs_limit: 220_000,
        },
        GpuTier::IntegratedOrLaptop => GpuAutoTuning {
            stream_matrix_budget_bytes: 128 * 1024 * 1024,
            max_obs_upload: 16_000_000,
            ready_batch_obs_limit: 120_000,
        },
    }
}

pub(crate) fn effective_matrix_budget(limits: &wgpu::Limits, requested_budget: usize) -> usize {
    requested_budget
        .min(limits.max_storage_buffer_binding_size as usize)
        .min(limits.max_buffer_size as usize)
        .max(1)
}

/// Cap max_obs_upload so the obs buffer never exceeds GPU binding/buffer-size limits and
/// a single dispatch_workgroups call never exceeds max_compute_workgroups_per_dimension.
/// obs_size_bytes  = size_of::<Observation>() — passed in to avoid importing that type here.
/// workgroup_size  = WORKGROUP_SIZE from kernel.rs.
pub(crate) fn effective_max_obs_upload(
    limits: &wgpu::Limits,
    requested: usize,
    obs_size_bytes: usize,
    workgroup_size: usize,
) -> usize {
    let max_from_binding =
        (limits.max_storage_buffer_binding_size as usize).min(limits.max_buffer_size as usize);
    let max_from_workgroups =
        (limits.max_compute_workgroups_per_dimension as usize).saturating_mul(workgroup_size);
    requested
        .min(max_from_binding / obs_size_bytes.max(1))
        .min(max_from_workgroups)
        .max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn adapter_info(name: &str, device_type: wgpu::DeviceType) -> wgpu::AdapterInfo {
        wgpu::AdapterInfo {
            name: name.to_string(),
            vendor: 0,
            device: 0,
            device_type,
            device_pci_bus_id: String::new(),
            driver: String::new(),
            driver_info: String::new(),
            backend: wgpu::Backend::Vulkan,
            subgroup_min_size: 0,
            subgroup_max_size: 0,
            transient_saves_memory: false,
        }
    }

    fn summary(index: usize, name: &str) -> GpuAdapterSummary {
        GpuAdapterSummary {
            index,
            info: adapter_info(name, wgpu::DeviceType::DiscreteGpu),
            limits: wgpu::Limits::default(),
        }
    }

    #[test]
    fn classify_adapter_detects_known_tiers() {
        assert_eq!(
            classify_adapter(&adapter_info(
                "NVIDIA H100 PCIe",
                wgpu::DeviceType::DiscreteGpu
            )),
            GpuTier::Datacenter
        );
        assert_eq!(
            classify_adapter(&adapter_info(
                "NVIDIA GeForce RTX 4090",
                wgpu::DeviceType::DiscreteGpu
            )),
            GpuTier::HighEndDiscrete
        );
        assert_eq!(
            classify_adapter(&adapter_info(
                "AMD Radeon Pro",
                wgpu::DeviceType::DiscreteGpu
            )),
            GpuTier::Discrete
        );
        assert_eq!(
            classify_adapter(&adapter_info(
                "Intel Iris Xe",
                wgpu::DeviceType::IntegratedGpu
            )),
            GpuTier::IntegratedOrLaptop
        );
    }

    #[test]
    fn auto_tuning_matches_expected_tier_defaults() {
        assert_eq!(
            auto_tuning_for_tier(GpuTier::Datacenter),
            GpuAutoTuning {
                stream_matrix_budget_bytes: 1024 * 1024 * 1024,
                max_obs_upload: 64_000_000,
                ready_batch_obs_limit: 750_000,
            }
        );
        assert_eq!(
            auto_tuning_for_tier(GpuTier::IntegratedOrLaptop),
            GpuAutoTuning {
                stream_matrix_budget_bytes: 128 * 1024 * 1024,
                max_obs_upload: 16_000_000,
                ready_batch_obs_limit: 120_000,
            }
        );
    }

    #[test]
    fn effective_matrix_budget_caps_to_wgpu_limits() {
        let limits = wgpu::Limits {
            max_storage_buffer_binding_size: 64 * 1024 * 1024,
            max_buffer_size: 128 * 1024 * 1024,
            ..Default::default()
        };

        assert_eq!(
            effective_matrix_budget(&limits, 256 * 1024 * 1024),
            64 * 1024 * 1024
        );
        assert_eq!(
            effective_matrix_budget(&limits, 32 * 1024 * 1024),
            32 * 1024 * 1024
        );
        assert_eq!(effective_matrix_budget(&limits, 0), 1);
    }

    #[test]
    fn effective_max_obs_upload_caps_by_binding_and_buffer_limits() {
        let limits = wgpu::Limits {
            max_storage_buffer_binding_size: 128 * 1024 * 1024, // 128 MiB
            max_buffer_size: 256 * 1024 * 1024,
            max_compute_workgroups_per_dimension: 65535,
            ..Default::default()
        };
        // 24M * 8 = 192 MiB > 128 MiB binding → capped to 16,777,216
        // but 16,777,216 / 256 = 65536 > 65535 workgroup limit → further capped to 65535*256 = 16,776,960
        assert_eq!(
            effective_max_obs_upload(&limits, 24_000_000, 8, 256),
            65535 * 256
        );
        // Already within limits: unchanged
        assert_eq!(
            effective_max_obs_upload(&limits, 1_000_000, 8, 256),
            1_000_000
        );
        // Zero obs_size_bytes: treated as 1-byte obs → workgroup limit dominates
        assert_eq!(
            effective_max_obs_upload(&limits, usize::MAX, 0, 256),
            65535 * 256
        );
    }

    #[test]
    fn select_adapter_index_uses_explicit_index() -> Result<()> {
        let summaries = vec![summary(0, "GPU 0"), summary(1, "GPU 1")];
        let selector = GpuSelector {
            index: Some(1),
            ..Default::default()
        };

        assert_eq!(select_adapter_index(&summaries, &selector)?, 1);
        Ok(())
    }

    #[test]
    fn select_adapter_index_uses_case_insensitive_name_match() -> Result<()> {
        let summaries = vec![summary(0, "Intel UHD"), summary(1, "NVIDIA H100 PCIe")];
        let selector = GpuSelector {
            name: Some("h100".to_string()),
            ..Default::default()
        };

        assert_eq!(select_adapter_index(&summaries, &selector)?, 1);
        Ok(())
    }

    #[test]
    fn select_adapter_index_rejects_conflicting_selectors() {
        let summaries = vec![summary(0, "GPU 0")];
        let selector = GpuSelector {
            index: Some(0),
            name: Some("gpu".to_string()),
            ..Default::default()
        };

        let err = select_adapter_index(&summaries, &selector).unwrap_err();
        assert!(err.to_string().contains("cannot be used together"));
    }

    #[test]
    fn select_adapter_index_rejects_out_of_range_index() {
        let summaries = vec![summary(0, "GPU 0")];
        let selector = GpuSelector {
            index: Some(2),
            ..Default::default()
        };

        let err = select_adapter_index(&summaries, &selector).unwrap_err();
        assert!(err.to_string().contains("out of range"));
    }
}
