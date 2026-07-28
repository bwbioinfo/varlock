use std::{collections::BTreeMap, fs::File, path::Path, time::Instant};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_bam::io::Reader;
use noodles_sam::alignment::{record::cigar::Op, record::cigar::op::Kind};

use super::samples::InputSampleResolver;
use super::targets::in_targets;
use super::types::{SiteCounts, SiteKey, TargetIndex, base_index};

pub(crate) struct PileupSettings<'a> {
    pub(crate) targets: &'a TargetIndex,
    pub(crate) sample_count: usize,
    pub(crate) min_baseq: u8,
}

pub(crate) struct ScanParams<'a> {
    pub(crate) sample_resolver: &'a InputSampleResolver,
    pub(crate) min_mapq: u8,
    pub(crate) verbose: u8,
    pub(crate) pileup: PileupSettings<'a>,
}

pub(crate) struct InputResult {
    pub(crate) counts: BTreeMap<SiteKey, SiteCounts>,
    pub(crate) skipped_rg: usize,
    pub(crate) skipped_flags: usize,
}

pub(crate) fn process_input_bam(path: &Path, scan: &ScanParams<'_>) -> Result<InputResult> {
    let started = Instant::now();
    if scan.verbose > 1 {
        eprintln!("[call_targets] scanning input {}", path.display());
    }

    let file =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = Reader::new(file);
    let _bam_header = reader
        .read_header()
        .with_context(|| format!("failed to read header for {}", path.display()))?;

    let mut counts: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
    let mut skipped_rg = 0usize;
    let mut skipped_flags = 0usize;
    let mut records_seen = 0u64;
    let mut records_retained = 0u64;
    let mut base_observations = 0u64;

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

        if !scan
            .pileup
            .targets
            .by_ref
            .contains_key(&reference_sequence_id)
        {
            continue;
        }

        let Some(sample_index) = scan.sample_resolver.resolve_record(&record)? else {
            skipped_rg += 1;
            continue;
        };

        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < scan.min_mapq {
            continue;
        }

        records_retained += 1;

        let added = pileup_record(
            &record,
            reference_sequence_id,
            sample_index,
            &mut counts,
            &scan.pileup,
        )?;
        base_observations += added as u64;
    }

    if scan.verbose > 1 {
        eprintln!(
            "[call_targets] finished input {} elapsed={:.2?} seen={} retained={} base_obs={} skipped_flags={} skipped_rg={}",
            path.display(),
            started.elapsed(),
            records_seen,
            records_retained,
            base_observations,
            skipped_flags,
            skipped_rg
        );
    }

    Ok(InputResult {
        counts,
        skipped_rg,
        skipped_flags,
    })
}

pub(crate) fn merge_counts(
    dst: &mut BTreeMap<SiteKey, SiteCounts>,
    src: BTreeMap<SiteKey, SiteCounts>,
) -> Result<()> {
    for (key, site_counts) in src {
        let sample_count = site_counts.per_sample.len();
        let dst_entry = dst.entry(key).or_insert_with(|| SiteCounts {
            per_sample: vec![[0; 4]; sample_count],
        });

        if dst_entry.per_sample.len() != sample_count {
            bail!("inconsistent sample vector length while merging counts");
        }

        for (dst_sample, src_sample) in dst_entry.per_sample.iter_mut().zip(site_counts.per_sample)
        {
            for i in 0..4 {
                dst_sample[i] = dst_sample[i].saturating_add(src_sample[i]);
            }
        }
    }

    Ok(())
}

#[cfg_attr(not(feature = "wgpu"), allow(dead_code))]
pub(crate) fn cap_counts(counts: &mut BTreeMap<SiteKey, SiteCounts>, max_depth: u32) {
    for site_counts in counts.values_mut() {
        for sample_counts in &mut site_counts.per_sample {
            let raw = *sample_counts;
            *sample_counts = [0; 4];
            merge_sample_counts_with_cap(sample_counts, raw, max_depth);
        }
    }
}

pub(crate) fn merge_sample_counts_with_cap(dst: &mut [u32; 4], src: [u32; 4], max_depth: u32) {
    if max_depth == 0 {
        return;
    }

    let max_depth = u64::from(max_depth);
    let dst_total = dst.iter().map(|&count| u64::from(count)).sum::<u64>();
    if dst_total >= max_depth {
        return;
    }

    let space = max_depth - dst_total;
    let src_total = src.iter().map(|&count| u64::from(count)).sum::<u64>();
    if src_total == 0 {
        return;
    }

    if src_total <= space {
        for i in 0..4 {
            dst[i] = dst[i].saturating_add(src[i]);
        }
        return;
    }

    // Deterministic proportional merge when source exceeds remaining max_depth.
    let mut add = [0u32; 4];
    let mut used = 0u64;
    let mut remainders = [(0u64, 0usize); 4];
    for i in 0..4 {
        let weighted = u64::from(src[i]) * space;
        add[i] = u32::try_from(weighted / src_total)
            .expect("proportional depth allocation must fit in u32");
        used += u64::from(add[i]);
        remainders[i] = (weighted % src_total, i);
    }

    let mut remaining = space - used;
    remainders.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    for &(_, i) in &remainders {
        if remaining == 0 {
            break;
        }
        if add[i] < src[i] {
            add[i] += 1;
            remaining -= 1;
        }
    }

    for i in 0..4 {
        dst[i] = dst[i].saturating_add(add[i]);
    }
}

pub(crate) fn should_skip_record(record: &bam::Record) -> bool {
    let flags = record.flags();
    flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_supplementary()
        || flags.is_qc_fail()
        || flags.is_duplicate()
}

fn pileup_record(
    record: &bam::Record,
    reference_sequence_id: usize,
    sample_index: usize,
    counts: &mut BTreeMap<SiteKey, SiteCounts>,
    settings: &PileupSettings<'_>,
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

    let intervals = match settings.targets.by_ref.get(&reference_sequence_id) {
        Some(v) => v,
        None => return Ok(0),
    };

    let mut ref_pos = alignment_start.get() as u64;
    let mut read_pos = 0usize;
    let mut added = 0u32;

    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let pos0 = ref_pos - 1;
                    if in_targets(intervals, pos0) {
                        let base = seq_buf.get(read_pos).unwrap_or(b'N');
                        let q = qual.get(read_pos).copied().unwrap_or(0);
                        if q >= settings.min_baseq
                            && let Some(idx) = base_index(base)
                        {
                            let key = SiteKey {
                                reference_sequence_id,
                                position: ref_pos as u32,
                            };
                            let entry = counts.entry(key).or_insert_with(|| SiteCounts {
                                per_sample: vec![[0; 4]; settings.sample_count],
                            });
                            let sample_counts = &mut entry.per_sample[sample_index];
                            sample_counts[idx] = sample_counts[idx].saturating_add(1);
                            added = added.saturating_add(1);
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
    use super::super::types::{SiteCounts, SiteKey};
    use super::{cap_counts, merge_counts, merge_sample_counts_with_cap};
    use anyhow::Result;
    use std::collections::BTreeMap;

    fn site(ref_id: usize, pos: u32) -> SiteKey {
        SiteKey {
            reference_sequence_id: ref_id,
            position: pos,
        }
    }

    fn counts(per_sample: Vec<[u32; 4]>) -> SiteCounts {
        SiteCounts { per_sample }
    }

    #[test]
    fn merge_sample_counts_adds_without_truncation() {
        let mut dst = [10, 0, 0, 0];
        merge_sample_counts_with_cap(&mut dst, [0, 5, 5, 0], 30);
        assert_eq!(dst, [10, 5, 5, 0]);
    }

    #[test]
    fn merge_sample_counts_truncates_deterministically() {
        let mut dst = [0, 0, 0, 0];
        merge_sample_counts_with_cap(&mut dst, [1, 1, 1, 1], 3);
        assert_eq!(dst, [1, 1, 1, 0]); // largest-remainder rounding
    }

    #[test]
    fn merge_sample_counts_at_capacity_adds_nothing() {
        let mut dst = [30, 0, 0, 0];
        merge_sample_counts_with_cap(&mut dst, [0, 5, 5, 0], 30);
        assert_eq!(dst, [30, 0, 0, 0]);
    }

    #[test]
    fn merge_sample_counts_zero_max_depth_adds_nothing() {
        let mut dst = [1, 2, 3, 4];
        let original = dst;
        merge_sample_counts_with_cap(&mut dst, [1, 1, 1, 1], 0);
        assert_eq!(dst, original);
    }

    #[test]
    fn merge_sample_counts_empty_source_adds_nothing() {
        let mut dst = [5, 5, 0, 0];
        let original = dst;
        merge_sample_counts_with_cap(&mut dst, [0, 0, 0, 0], 100);
        assert_eq!(dst, original);
    }

    #[test]
    fn merge_counts_accumulates_across_calls() -> Result<()> {
        let mut all: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut first = BTreeMap::new();
        first.insert(site(0, 100), counts(vec![[5, 0, 0, 0]]));

        let mut second = BTreeMap::new();
        second.insert(site(0, 100), counts(vec![[3, 0, 0, 0]]));

        merge_counts(&mut all, first)?;
        merge_counts(&mut all, second)?;

        assert_eq!(all[&site(0, 100)].per_sample[0], [8, 0, 0, 0]);
        Ok(())
    }

    #[test]
    fn merge_counts_inserts_new_sites() -> Result<()> {
        let mut all: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut src = BTreeMap::new();
        src.insert(site(0, 1), counts(vec![[0, 5, 0, 0]]));
        src.insert(site(0, 2), counts(vec![[0, 0, 7, 0]]));

        merge_counts(&mut all, src)?;

        assert_eq!(all.len(), 2);
        assert_eq!(all[&site(0, 1)].per_sample[0][1], 5);
        assert_eq!(all[&site(0, 2)].per_sample[0][2], 7);
        Ok(())
    }

    #[test]
    fn merge_counts_errors_on_sample_count_mismatch() -> Result<()> {
        let mut all: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut first = BTreeMap::new();
        first.insert(site(0, 1), counts(vec![[1, 0, 0, 0], [0, 1, 0, 0]]));
        merge_counts(&mut all, first)?;

        let mut second = BTreeMap::new();
        second.insert(site(0, 1), counts(vec![[0, 0, 1, 0]])); // 1 sample vs 2
        let err = merge_counts(&mut all, second).unwrap_err();
        assert!(err.to_string().contains("inconsistent sample vector"));
        Ok(())
    }

    #[test]
    fn uncapped_merge_then_cap_matches_single_cap() -> Result<()> {
        let mut split: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut first = BTreeMap::new();
        first.insert(site(0, 1), counts(vec![[0, 5, 1, 0]]));
        merge_counts(&mut split, first)?;

        let mut second = BTreeMap::new();
        second.insert(site(0, 1), counts(vec![[0, 5, 3, 0]]));
        merge_counts(&mut split, second)?;
        cap_counts(&mut split, 6);

        let mut single: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();
        single.insert(site(0, 1), counts(vec![[0, 10, 4, 0]]));
        cap_counts(&mut single, 6);

        assert_eq!(
            split[&site(0, 1)].per_sample[0],
            single[&site(0, 1)].per_sample[0]
        );
        Ok(())
    }

    #[test]
    fn merge_counts_then_cap_is_input_order_independent() -> Result<()> {
        let first = counts(vec![[6, 0, 0, 0]]);
        let second = counts(vec![[0, 6, 0, 0]]);

        let mut forward = BTreeMap::new();
        forward.insert(site(0, 1), first.clone());
        merge_counts(&mut forward, BTreeMap::from([(site(0, 1), second.clone())]))?;
        cap_counts(&mut forward, 5);

        let mut reverse = BTreeMap::new();
        reverse.insert(site(0, 1), second);
        merge_counts(&mut reverse, BTreeMap::from([(site(0, 1), first)]))?;
        cap_counts(&mut reverse, 5);

        assert_eq!(
            forward[&site(0, 1)].per_sample,
            reverse[&site(0, 1)].per_sample
        );
        assert_eq!(forward[&site(0, 1)].per_sample[0], [3, 2, 0, 0]);
        Ok(())
    }

    #[test]
    fn merge_counts_saturates_raw_base_counters() -> Result<()> {
        let mut merged = BTreeMap::from([(site(0, 1), counts(vec![[u32::MAX, 0, 0, 0]]))]);
        merge_counts(
            &mut merged,
            BTreeMap::from([(site(0, 1), counts(vec![[1, 0, 0, 0]]))]),
        )?;

        assert_eq!(merged[&site(0, 1)].per_sample[0], [u32::MAX, 0, 0, 0]);
        Ok(())
    }
}
