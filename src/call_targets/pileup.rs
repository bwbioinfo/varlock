use std::{collections::BTreeMap, fs::File, path::Path, time::Instant};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_bam::io::Reader;
use noodles_sam::alignment::{record::cigar::Op, record::cigar::op::Kind};

use super::samples::InputSampleResolver;
use super::targets::in_targets;
use super::types::{
    IndelAllele, IndelCounts, IndelKey, SiteCounts, SiteKey, TargetIndex, base_index,
};

pub(crate) struct PileupSettings<'a> {
    pub(crate) targets: &'a TargetIndex,
    pub(crate) sample_count: usize,
    pub(crate) min_baseq: u8,
    pub(crate) max_depth: u32,
}

pub(crate) struct ScanParams<'a> {
    pub(crate) sample_resolver: &'a InputSampleResolver,
    pub(crate) min_mapq: u8,
    pub(crate) verbose: u8,
    pub(crate) pileup: PileupSettings<'a>,
}

pub(crate) struct InputResult {
    pub(crate) counts: BTreeMap<SiteKey, SiteCounts>,
    pub(crate) indel_counts: BTreeMap<IndelKey, IndelCounts>,
    pub(crate) skipped_rg: usize,
    pub(crate) skipped_flags: usize,
}

pub(crate) fn process_input_bam(
    path: &Path,
    scan: &ScanParams<'_>,
    collect_indels: bool,
) -> Result<InputResult> {
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
    let mut indel_alt_counts = collect_indels.then(BTreeMap::new);
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
            indel_alt_counts.as_mut(),
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
        indel_counts: finalize_indel_counts(indel_alt_counts.unwrap_or_default(), &counts)?,
        counts,
        skipped_rg,
        skipped_flags,
    })
}

/// Collects CIGAR-derived indel alternate support without retaining whole-genome base counts.
/// The GPU path uses this after device aggregation has produced anchor depth counts.
#[cfg(feature = "wgpu")]
pub(crate) fn collect_indel_alt_counts(
    path: &Path,
    scan: &ScanParams<'_>,
) -> Result<BTreeMap<IndelKey, Vec<u32>>> {
    let file =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = Reader::new(file);
    let _bam_header = reader
        .read_header()
        .with_context(|| format!("failed to read header for {}", path.display()))?;

    let mut counts = BTreeMap::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("failed to read record {}", path.display()))?;
        if should_skip_record(&record) {
            continue;
        }
        let reference_sequence_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            Some(Err(e)) => return Err(e).context("failed to read reference sequence id"),
            None => continue,
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
            continue;
        };
        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < scan.min_mapq {
            continue;
        }
        collect_record_indel_alts(
            &record,
            reference_sequence_id,
            sample_index,
            &scan.pileup,
            &mut counts,
        )?;
    }

    Ok(counts)
}

pub(crate) fn merge_counts(
    dst: &mut BTreeMap<SiteKey, SiteCounts>,
    src: BTreeMap<SiteKey, SiteCounts>,
    max_depth: u32,
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
            merge_sample_counts_with_cap(dst_sample, src_sample, max_depth);
        }
    }

    Ok(())
}

pub(crate) fn merge_indel_counts(
    dst: &mut BTreeMap<IndelKey, IndelCounts>,
    src: BTreeMap<IndelKey, IndelCounts>,
    max_depth: u32,
) -> Result<()> {
    for (key, indel_counts) in src {
        let sample_count = indel_counts.per_sample.len();
        let dst_entry = dst.entry(key).or_insert_with(|| IndelCounts {
            per_sample: vec![[0; 2]; sample_count],
        });
        if dst_entry.per_sample.len() != sample_count {
            bail!("inconsistent sample vector length while merging indel counts");
        }
        for (dst_sample, src_sample) in dst_entry.per_sample.iter_mut().zip(indel_counts.per_sample)
        {
            merge_allele_counts_with_cap(dst_sample, src_sample, max_depth);
        }
    }

    Ok(())
}

#[cfg(feature = "wgpu")]
pub(crate) fn merge_indel_alt_counts(
    dst: &mut BTreeMap<IndelKey, Vec<u32>>,
    src: BTreeMap<IndelKey, Vec<u32>>,
) -> Result<()> {
    for (key, src_counts) in src {
        let sample_count = src_counts.len();
        let dst_counts = dst.entry(key).or_insert_with(|| vec![0; sample_count]);
        if dst_counts.len() != sample_count {
            bail!("inconsistent sample vector length while merging indel alternate counts");
        }
        for (dst_count, src_count) in dst_counts.iter_mut().zip(src_counts) {
            *dst_count = dst_count.saturating_add(src_count);
        }
    }
    Ok(())
}

pub(crate) fn finalize_indel_counts(
    indel_alt_counts: BTreeMap<IndelKey, Vec<u32>>,
    site_counts: &BTreeMap<SiteKey, SiteCounts>,
) -> Result<BTreeMap<IndelKey, IndelCounts>> {
    let mut indel_counts = BTreeMap::new();
    for (key, alt_counts) in indel_alt_counts {
        let anchor_counts = site_counts
            .get(&key.anchor_site())
            .context("indel evidence is missing its anchor depth count")?;
        if anchor_counts.per_sample.len() != alt_counts.len() {
            bail!("inconsistent sample vector length while finalizing indel counts");
        }
        let per_sample = anchor_counts
            .per_sample
            .iter()
            .zip(alt_counts)
            .map(|(base_counts, alt_count)| {
                let depth = base_counts.iter().sum::<u32>();
                let alt_count = alt_count.min(depth);
                [depth - alt_count, alt_count]
            })
            .collect();
        indel_counts.insert(key, IndelCounts { per_sample });
    }
    Ok(indel_counts)
}

#[cfg_attr(not(feature = "wgpu"), allow(dead_code))]
pub(crate) fn merge_counts_uncapped(
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
    merge_allele_counts_with_cap(dst, src, max_depth);
}

fn merge_allele_counts_with_cap<const N: usize>(dst: &mut [u32; N], src: [u32; N], max_depth: u32) {
    if max_depth == 0 {
        return;
    }

    let dst_total = dst.iter().sum::<u32>();
    if dst_total >= max_depth {
        return;
    }

    let space = max_depth - dst_total;
    let src_total = src.iter().sum::<u32>();
    if src_total == 0 {
        return;
    }

    if src_total <= space {
        for i in 0..N {
            dst[i] += src[i];
        }
        return;
    }

    // Deterministic proportional merge when source exceeds remaining max_depth.
    let mut add = [0u32; N];
    let mut used = 0u32;
    let mut remainders = [(0u64, 0usize); N];
    for i in 0..N {
        let weighted = src[i] as u64 * space as u64;
        add[i] = (weighted / src_total as u64) as u32;
        used += add[i];
        remainders[i] = (weighted % src_total as u64, i);
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

    for i in 0..N {
        dst[i] += add[i];
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
    mut indel_alt_counts: Option<&mut BTreeMap<IndelKey, Vec<u32>>>,
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
    let sequence = (0..seq_buf.len())
        .map(|index| seq_buf.get(index).unwrap_or(b'N'))
        .collect::<Vec<_>>();

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
    let mut last_anchor = None;

    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let pos0 = ref_pos - 1;
                    let mut accepted_anchor = None;
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
                            let dp = sample_counts.iter().sum::<u32>();
                            if dp < settings.max_depth {
                                sample_counts[idx] += 1;
                                added += 1;
                                accepted_anchor = u32::try_from(ref_pos).ok();
                            }
                        }
                    }
                    last_anchor = accepted_anchor;
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion => {
                if let Some(anchor_position) =
                    last_anchor.filter(|&position| u64::from(position) == ref_pos.saturating_sub(1))
                    && let Some(indel_alt_counts) = indel_alt_counts.as_deref_mut()
                    && let Some(inserted) =
                        insertion_allele(&sequence, qual, read_pos, len, settings.min_baseq)
                {
                    insert_indel_alt(
                        indel_alt_counts,
                        IndelKey {
                            reference_sequence_id,
                            position: anchor_position,
                            allele: IndelAllele::Insertion(inserted),
                        },
                        sample_index,
                        settings.sample_count,
                    )?;
                }
                read_pos += len;
            }
            Kind::SoftClip => {
                read_pos += len;
                last_anchor = None;
            }
            Kind::Deletion => {
                if let Some(anchor_position) =
                    last_anchor.filter(|&position| u64::from(position) == ref_pos.saturating_sub(1))
                    && let Some(indel_alt_counts) = indel_alt_counts.as_deref_mut()
                {
                    let deletion =
                        u32::try_from(len).context("CIGAR deletion length exceeds u32")?;
                    insert_indel_alt(
                        indel_alt_counts,
                        IndelKey {
                            reference_sequence_id,
                            position: anchor_position,
                            allele: IndelAllele::Deletion(deletion),
                        },
                        sample_index,
                        settings.sample_count,
                    )?;
                }
                ref_pos += len as u64;
                last_anchor = None;
            }
            Kind::Skip => {
                ref_pos += len as u64;
                last_anchor = None;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(added)
}

#[cfg(feature = "wgpu")]
fn collect_record_indel_alts(
    record: &bam::Record,
    reference_sequence_id: usize,
    sample_index: usize,
    settings: &PileupSettings<'_>,
    indel_alt_counts: &mut BTreeMap<IndelKey, Vec<u32>>,
) -> Result<()> {
    let alignment_start = match record.alignment_start() {
        Some(Ok(pos)) => pos,
        Some(Err(e)) => return Err(e).context("failed to read alignment start"),
        None => return Ok(()),
    };
    let seq_buf = record.sequence();
    let qual_buf = record.quality_scores();
    let qual = qual_buf.as_ref();
    if seq_buf.len() != qual.len() {
        return Ok(());
    }
    let sequence = (0..seq_buf.len())
        .map(|index| seq_buf.get(index).unwrap_or(b'N'))
        .collect::<Vec<_>>();
    let ops: Vec<Op> = record
        .cigar()
        .iter()
        .collect::<std::io::Result<Vec<_>>>()
        .context("failed to read CIGAR")?;
    let intervals = match settings.targets.by_ref.get(&reference_sequence_id) {
        Some(v) => v,
        None => return Ok(()),
    };

    let mut ref_pos = alignment_start.get() as u64;
    let mut read_pos = 0usize;
    let mut last_anchor = None;
    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let base = seq_buf.get(read_pos).unwrap_or(b'N');
                    let q = qual.get(read_pos).copied().unwrap_or(0);
                    last_anchor = (in_targets(intervals, ref_pos - 1)
                        && q >= settings.min_baseq
                        && base_index(base).is_some())
                    .then(|| u32::try_from(ref_pos).ok())
                    .flatten();
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion => {
                if let Some(anchor_position) =
                    last_anchor.filter(|&position| u64::from(position) == ref_pos.saturating_sub(1))
                    && let Some(inserted) =
                        insertion_allele(&sequence, qual, read_pos, len, settings.min_baseq)
                {
                    insert_indel_alt(
                        indel_alt_counts,
                        IndelKey {
                            reference_sequence_id,
                            position: anchor_position,
                            allele: IndelAllele::Insertion(inserted),
                        },
                        sample_index,
                        settings.sample_count,
                    )?;
                }
                read_pos += len;
            }
            Kind::SoftClip => {
                read_pos += len;
                last_anchor = None;
            }
            Kind::Deletion => {
                if let Some(anchor_position) =
                    last_anchor.filter(|&position| u64::from(position) == ref_pos.saturating_sub(1))
                {
                    let deletion =
                        u32::try_from(len).context("CIGAR deletion length exceeds u32")?;
                    insert_indel_alt(
                        indel_alt_counts,
                        IndelKey {
                            reference_sequence_id,
                            position: anchor_position,
                            allele: IndelAllele::Deletion(deletion),
                        },
                        sample_index,
                        settings.sample_count,
                    )?;
                }
                ref_pos += len as u64;
                last_anchor = None;
            }
            Kind::Skip => {
                ref_pos += len as u64;
                last_anchor = None;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(())
}

fn insertion_allele(
    sequence: &[u8],
    qualities: &[u8],
    read_pos: usize,
    len: usize,
    min_baseq: u8,
) -> Option<Vec<u8>> {
    let end = read_pos.checked_add(len)?;
    let (Some(bases), Some(quality_scores)) =
        (sequence.get(read_pos..end), qualities.get(read_pos..end))
    else {
        return None;
    };
    let inserted = bases
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    (quality_scores.iter().all(|&q| q >= min_baseq)
        && inserted.iter().all(|&base| base_index(base).is_some()))
    .then_some(inserted)
}

fn insert_indel_alt(
    indel_alt_counts: &mut BTreeMap<IndelKey, Vec<u32>>,
    key: IndelKey,
    sample_index: usize,
    sample_count: usize,
) -> Result<()> {
    let counts = indel_alt_counts
        .entry(key)
        .or_insert_with(|| vec![0; sample_count]);
    let count = counts
        .get_mut(sample_index)
        .context("indel sample index out of range")?;
    *count = count.saturating_add(1);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::types::{SiteCounts, SiteKey};
    use super::{cap_counts, merge_counts, merge_counts_uncapped, merge_sample_counts_with_cap};
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

        merge_counts(&mut all, first, 1000)?;
        merge_counts(&mut all, second, 1000)?;

        assert_eq!(all[&site(0, 100)].per_sample[0], [8, 0, 0, 0]);
        Ok(())
    }

    #[test]
    fn merge_counts_inserts_new_sites() -> Result<()> {
        let mut all: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut src = BTreeMap::new();
        src.insert(site(0, 1), counts(vec![[0, 5, 0, 0]]));
        src.insert(site(0, 2), counts(vec![[0, 0, 7, 0]]));

        merge_counts(&mut all, src, 1000)?;

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
        merge_counts(&mut all, first, 1000)?;

        let mut second = BTreeMap::new();
        second.insert(site(0, 1), counts(vec![[0, 0, 1, 0]])); // 1 sample vs 2
        let err = merge_counts(&mut all, second, 1000).unwrap_err();
        assert!(err.to_string().contains("inconsistent sample vector"));
        Ok(())
    }

    #[test]
    fn uncapped_merge_then_cap_matches_single_cap() -> Result<()> {
        let mut split: BTreeMap<SiteKey, SiteCounts> = BTreeMap::new();

        let mut first = BTreeMap::new();
        first.insert(site(0, 1), counts(vec![[0, 5, 1, 0]]));
        merge_counts_uncapped(&mut split, first)?;

        let mut second = BTreeMap::new();
        second.insert(site(0, 1), counts(vec![[0, 5, 3, 0]]));
        merge_counts_uncapped(&mut split, second)?;
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
}
