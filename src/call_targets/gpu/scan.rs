#![allow(dead_code)]

use std::{
    fs::File,
    path::{Path, PathBuf},
    sync::{
        Arc,
        mpsc::{Receiver, SyncSender, TrySendError, sync_channel},
    },
    thread,
    time::Instant,
};

use anyhow::{Context, Result, anyhow};
use noodles_bam as bam;
use noodles_bam::io::Reader;
use noodles_sam::alignment::{record::cigar::Op, record::cigar::op::Kind};

use crate::call_targets::{
    observation::{Observation, TargetFrontierIndex, TargetSiteMap},
    pileup::should_skip_record,
    samples::InputSampleResolver,
    types::{SiteKey, base_index},
};

const SCAN_BATCH_OBSERVATIONS: usize = 100_000;
const FRONTIER_PROGRESS_STEP: u32 = 20_000;

#[derive(Debug)]
pub(crate) enum ScanEvent {
    Batch {
        input_idx: usize,
        frontier_site_idx: u32,
        observations: Vec<Observation>,
    },
    Progress {
        input_idx: usize,
        frontier_site_idx: u32,
    },
    Done {
        input_idx: usize,
        skipped_rg: usize,
        skipped_flags: usize,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CoveredObservation {
    pub(crate) site_key: SiteKey,
    pub(crate) sample_base: u32,
}

#[derive(Debug)]
pub(crate) enum CoveredScanEvent {
    Batch {
        input_idx: usize,
        observations: Vec<CoveredObservation>,
    },
    Done {
        input_idx: usize,
        skipped_rg: usize,
        skipped_flags: usize,
    },
}

#[derive(Clone)]
pub(crate) struct ScanWorkerParams {
    pub(crate) input_sample_resolvers: Arc<Vec<InputSampleResolver>>,
    pub(crate) min_mapq: u8,
    pub(crate) min_baseq: u8,
    pub(crate) max_depth: u32,
    pub(crate) verbose: u8,
}

pub(crate) fn spawn_scan_workers(
    inputs: &[PathBuf],
    params: ScanWorkerParams,
    target_sites: Arc<TargetSiteMap>,
    frontier_index: Arc<TargetFrontierIndex>,
) -> (Receiver<ScanEvent>, Vec<thread::JoinHandle<Result<()>>>) {
    let capacity = scan_channel_capacity(inputs.len());
    let (scan_tx, scan_rx) = sync_channel(capacity);
    let params = Arc::new(params);

    let handles = inputs
        .iter()
        .cloned()
        .enumerate()
        .map(|(input_idx, path)| {
            let params = Arc::clone(&params);
            let target_sites = Arc::clone(&target_sites);
            let frontier_index = Arc::clone(&frontier_index);
            let scan_tx = scan_tx.clone();
            thread::spawn(move || {
                scan_bam_worker(
                    input_idx,
                    &path,
                    &params,
                    &target_sites,
                    &frontier_index,
                    scan_tx,
                )
            })
        })
        .collect();

    (scan_rx, handles)
}

pub(crate) fn spawn_covered_scan_workers(
    inputs: &[PathBuf],
    params: ScanWorkerParams,
) -> (
    Receiver<CoveredScanEvent>,
    Vec<thread::JoinHandle<Result<()>>>,
) {
    let capacity = scan_channel_capacity(inputs.len());
    let (scan_tx, scan_rx) = sync_channel(capacity);
    let params = Arc::new(params);

    let handles = inputs
        .iter()
        .cloned()
        .enumerate()
        .map(|(input_idx, path)| {
            let params = Arc::clone(&params);
            let scan_tx = scan_tx.clone();
            thread::spawn(move || scan_bam_covered_worker(input_idx, &path, &params, scan_tx))
        })
        .collect();

    (scan_rx, handles)
}

fn scan_channel_capacity(input_count: usize) -> usize {
    32.max(input_count.saturating_mul(4))
}

pub(crate) fn scan_bam_worker(
    input_idx: usize,
    path: &Path,
    params: &ScanWorkerParams,
    target_sites: &TargetSiteMap,
    frontier_index: &TargetFrontierIndex,
    scan_tx: SyncSender<ScanEvent>,
) -> Result<()> {
    let sample_resolver = params
        .input_sample_resolvers
        .get(input_idx)
        .with_context(|| format!("missing sample resolver for input {}", path.display()))?;
    let started = Instant::now();
    if params.verbose > 1 {
        eprintln!("[call_targets_gpu] scanning input {}", path.display());
    }

    let file =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = Reader::new(file);
    let _bam_header = reader
        .read_header()
        .with_context(|| format!("failed to read header for {}", path.display()))?;

    let mut observations = Vec::with_capacity(SCAN_BATCH_OBSERVATIONS);
    let mut skipped_rg = 0usize;
    let mut skipped_flags = 0usize;
    let mut records_seen = 0u64;
    let mut records_retained = 0u64;
    let mut observation_count = 0u64;
    let mut frontier_site_idx = 0u32;
    let mut last_progress_frontier = 0u32;

    for result in reader.records() {
        let record = result.with_context(|| format!("failed to read record {}", path.display()))?;
        records_seen += 1;

        if should_skip_record(&record) {
            skipped_flags += 1;
            continue;
        }

        let reference_sequence_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            Some(Err(e)) => return Err(e).context("failed to read reference sequence id"),
            None => {
                skipped_flags += 1;
                continue;
            }
        };
        if let Some(Ok(position)) = record.alignment_start() {
            let position = u32::try_from(position.get()).unwrap_or(u32::MAX);
            let candidate = frontier_index.frontier_for_record(
                &target_sites.site_keys,
                reference_sequence_id,
                position,
            );
            frontier_site_idx = frontier_site_idx.max(candidate);
            if frontier_site_idx.saturating_sub(last_progress_frontier) >= FRONTIER_PROGRESS_STEP {
                send_scan_event(
                    &scan_tx,
                    ScanEvent::Progress {
                        input_idx,
                        frontier_site_idx,
                    },
                )?;
                last_progress_frontier = frontier_site_idx;
            }
        }

        if !target_sites.contains_ref(reference_sequence_id) {
            continue;
        }

        let Some(sample_index) = sample_resolver.resolve_record(&record)? else {
            skipped_rg += 1;
            continue;
        };

        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < params.min_mapq {
            continue;
        }

        records_retained += 1;
        let added = pileup_record_observations(
            &record,
            reference_sequence_id,
            sample_index,
            target_sites,
            params.min_baseq,
            &mut observations,
        )?;
        observation_count += added as u64;

        if observations.len() >= SCAN_BATCH_OBSERVATIONS {
            send_scan_event(
                &scan_tx,
                ScanEvent::Batch {
                    input_idx,
                    frontier_site_idx,
                    observations: std::mem::take(&mut observations),
                },
            )?;
        }
    }

    if !observations.is_empty() {
        send_scan_event(
            &scan_tx,
            ScanEvent::Batch {
                input_idx,
                frontier_site_idx,
                observations,
            },
        )?;
    } else if frontier_site_idx > last_progress_frontier {
        send_scan_event(
            &scan_tx,
            ScanEvent::Progress {
                input_idx,
                frontier_site_idx,
            },
        )?;
    }

    if params.verbose > 1 {
        eprintln!(
            "[call_targets_gpu] finished input {} elapsed={:.2?} seen={} retained={} observations={} skipped_flags={} skipped_rg={} max_depth={}",
            path.display(),
            started.elapsed(),
            records_seen,
            records_retained,
            observation_count,
            skipped_flags,
            skipped_rg,
            params.max_depth
        );
    }
    send_scan_event(
        &scan_tx,
        ScanEvent::Done {
            input_idx,
            skipped_rg,
            skipped_flags,
        },
    )?;

    Ok(())
}

pub(crate) fn scan_bam_covered_worker(
    input_idx: usize,
    path: &Path,
    params: &ScanWorkerParams,
    scan_tx: SyncSender<CoveredScanEvent>,
) -> Result<()> {
    let sample_resolver = params
        .input_sample_resolvers
        .get(input_idx)
        .with_context(|| format!("missing sample resolver for input {}", path.display()))?;
    let started = Instant::now();
    if params.verbose > 1 {
        eprintln!("[call_targets_gpu] scanning input {}", path.display());
    }

    let file =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = Reader::new(file);
    let _bam_header = reader
        .read_header()
        .with_context(|| format!("failed to read header for {}", path.display()))?;

    let mut observations = Vec::with_capacity(SCAN_BATCH_OBSERVATIONS);
    let mut skipped_rg = 0usize;
    let mut skipped_flags = 0usize;
    let mut records_seen = 0u64;
    let mut records_retained = 0u64;
    let mut observation_count = 0u64;

    for result in reader.records() {
        let record = result.with_context(|| format!("failed to read record {}", path.display()))?;
        records_seen += 1;

        if should_skip_record(&record) {
            skipped_flags += 1;
            continue;
        }

        let reference_sequence_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            Some(Err(e)) => return Err(e).context("failed to read reference sequence id"),
            None => {
                skipped_flags += 1;
                continue;
            }
        };

        let Some(sample_index) = sample_resolver.resolve_record(&record)? else {
            skipped_rg += 1;
            continue;
        };

        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < params.min_mapq {
            continue;
        }

        records_retained += 1;
        let added = pileup_record_covered_observations(
            &record,
            reference_sequence_id,
            sample_index,
            params.min_baseq,
            &mut observations,
        )?;
        observation_count += added as u64;

        if observations.len() >= SCAN_BATCH_OBSERVATIONS {
            send_covered_scan_event(
                &scan_tx,
                CoveredScanEvent::Batch {
                    input_idx,
                    observations: std::mem::take(&mut observations),
                },
            )?;
        }
    }

    if !observations.is_empty() {
        send_covered_scan_event(
            &scan_tx,
            CoveredScanEvent::Batch {
                input_idx,
                observations,
            },
        )?;
    }

    if params.verbose > 1 {
        eprintln!(
            "[call_targets_gpu] finished input {} elapsed={:.2?} seen={} retained={} observations={} skipped_flags={} skipped_rg={} max_depth={}",
            path.display(),
            started.elapsed(),
            records_seen,
            records_retained,
            observation_count,
            skipped_flags,
            skipped_rg,
            params.max_depth
        );
    }
    send_covered_scan_event(
        &scan_tx,
        CoveredScanEvent::Done {
            input_idx,
            skipped_rg,
            skipped_flags,
        },
    )?;

    Ok(())
}

fn send_scan_event(scan_tx: &SyncSender<ScanEvent>, event: ScanEvent) -> Result<()> {
    if matches!(event, ScanEvent::Progress { .. }) {
        match scan_tx.try_send(event) {
            Ok(()) | Err(TrySendError::Full(_)) => return Ok(()),
            Err(TrySendError::Disconnected(_)) => {
                return Err(anyhow!(
                    "call_targets_gpu aggregate worker stopped receiving scan events"
                ));
            }
        }
    }

    scan_tx
        .send(event)
        .map_err(|_| anyhow!("call_targets_gpu aggregate worker stopped receiving scan events"))
}

fn send_covered_scan_event(
    scan_tx: &SyncSender<CoveredScanEvent>,
    event: CoveredScanEvent,
) -> Result<()> {
    scan_tx
        .send(event)
        .map_err(|_| anyhow!("call_targets_gpu aggregate worker stopped receiving scan events"))
}

pub(crate) fn pileup_record_observations(
    record: &bam::Record,
    reference_sequence_id: usize,
    sample_index: usize,
    target_sites: &TargetSiteMap,
    min_baseq: u8,
    observations: &mut Vec<Observation>,
) -> Result<u32> {
    let alignment_start = match record.alignment_start() {
        Some(Ok(pos)) => pos,
        Some(Err(e)) => return Err(e).context("failed to read alignment start"),
        None => return Ok(0),
    };

    let seq_buf = record.sequence();
    let qual_buf = record.quality_scores();
    let qual = qual_buf.as_ref();
    if seq_buf.len() != qual.len() {
        return Ok(0);
    }

    let ops: Vec<Op> = record
        .cigar()
        .iter()
        .collect::<std::io::Result<Vec<_>>>()
        .context("failed to read CIGAR")?;

    let mut ref_pos = alignment_start.get() as u64;
    let mut read_pos = 0usize;
    let mut added = 0u32;

    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if ref_pos <= u32::MAX as u64 {
                        let pos1 = ref_pos as u32;
                        if let Some(site_idx) = target_sites.site_index(reference_sequence_id, pos1)
                        {
                            let base = seq_buf.get(read_pos).unwrap_or(b'N');
                            let q = qual.get(read_pos).copied().unwrap_or(0);
                            if q >= min_baseq
                                && let Some(base_idx) = base_index(base)
                            {
                                observations.push(Observation::new(
                                    site_idx,
                                    sample_index as u32,
                                    base_idx as u32,
                                )?);
                                added += 1;
                            }
                        }
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(added)
}

pub(crate) fn pileup_record_covered_observations(
    record: &bam::Record,
    reference_sequence_id: usize,
    sample_index: usize,
    min_baseq: u8,
    observations: &mut Vec<CoveredObservation>,
) -> Result<u32> {
    let alignment_start = match record.alignment_start() {
        Some(Ok(pos)) => pos,
        Some(Err(e)) => return Err(e).context("failed to read alignment start"),
        None => return Ok(0),
    };

    let seq_buf = record.sequence();
    let qual_buf = record.quality_scores();
    let qual = qual_buf.as_ref();
    if seq_buf.len() != qual.len() {
        return Ok(0);
    }

    let ops: Vec<Op> = record
        .cigar()
        .iter()
        .collect::<std::io::Result<Vec<_>>>()
        .context("failed to read CIGAR")?;

    let mut ref_pos = alignment_start.get() as u64;
    let mut read_pos = 0usize;
    let mut added = 0u32;

    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if ref_pos <= u32::MAX as u64 {
                        let base = seq_buf.get(read_pos).unwrap_or(b'N');
                        let q = qual.get(read_pos).copied().unwrap_or(0);
                        if q >= min_baseq
                            && let Some(base_idx) = base_index(base)
                        {
                            observations.push(CoveredObservation {
                                site_key: SiteKey {
                                    reference_sequence_id,
                                    position: ref_pos as u32,
                                },
                                sample_base: ((sample_index as u32) << 2) | base_idx as u32,
                            });
                            added += 1;
                        }
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(added)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_channel_capacity_has_floor_and_scales_with_inputs() {
        assert_eq!(scan_channel_capacity(0), 32);
        assert_eq!(scan_channel_capacity(4), 32);
        assert_eq!(scan_channel_capacity(9), 36);
    }

    #[test]
    fn send_scan_event_drops_progress_when_channel_is_full() -> Result<()> {
        let (tx, _rx) = sync_channel(0);

        send_scan_event(
            &tx,
            ScanEvent::Progress {
                input_idx: 0,
                frontier_site_idx: 10,
            },
        )?;

        Ok(())
    }
}
