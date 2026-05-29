use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::{Context, Result};

use super::types::{Interval, TargetIndex};

pub(crate) fn load_targets(
    path: &Path,
    ref_name_to_id: &HashMap<String, usize>,
) -> Result<TargetIndex> {
    let file = File::open(path)
        .with_context(|| format!("failed to open targets BED {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut by_ref: HashMap<usize, Vec<Interval>> = HashMap::new();

    for line in reader.lines() {
        let line = line.context("failed to read BED line")?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields.next().context("missing BED chrom")?;
        let start = fields
            .next()
            .context("missing BED start")?
            .parse::<u64>()
            .context("invalid BED start")?;
        let end = fields
            .next()
            .context("missing BED end")?
            .parse::<u64>()
            .context("invalid BED end")?;
        if end <= start {
            continue;
        }

        let reference_sequence_id = ref_name_to_id
            .get(chrom)
            .copied()
            .with_context(|| format!("BED chrom {} not in BAM header", chrom))?;
        by_ref
            .entry(reference_sequence_id)
            .or_default()
            .push(Interval { start, end });
    }

    for intervals in by_ref.values_mut() {
        intervals.sort_by_key(|i| i.start);
        let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
        for interval in intervals.drain(..) {
            if let Some(last) = merged.last_mut()
                && interval.start <= last.end
            {
                last.end = last.end.max(interval.end);
                continue;
            }
            merged.push(interval);
        }
        *intervals = merged;
    }

    Ok(TargetIndex { by_ref })
}

pub(crate) fn in_targets(intervals: &[Interval], pos0: u64) -> bool {
    let mut lo = 0usize;
    let mut hi = intervals.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        let interval = &intervals[mid];
        if pos0 < interval.start {
            hi = mid;
        } else if pos0 >= interval.end {
            lo = mid + 1;
        } else {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, io::Write};
    use anyhow::Result;
    use tempfile::NamedTempFile;
    use super::{in_targets, load_targets};
    use super::super::types::Interval;

    #[test]
    fn in_targets_hit_inside_interval() {
        let intervals = vec![Interval { start: 10, end: 20 }];
        assert!(in_targets(&intervals, 15));
    }

    #[test]
    fn in_targets_start_is_inclusive() {
        let intervals = vec![Interval { start: 10, end: 20 }];
        assert!(in_targets(&intervals, 10));
    }

    #[test]
    fn in_targets_end_is_exclusive() {
        let intervals = vec![Interval { start: 10, end: 20 }];
        assert!(!in_targets(&intervals, 20));
    }

    #[test]
    fn in_targets_misses_gap_between_intervals() {
        let intervals = vec![
            Interval { start: 10, end: 20 },
            Interval { start: 30, end: 40 },
        ];
        assert!(!in_targets(&intervals, 25));
    }

    #[test]
    fn in_targets_hits_second_interval() {
        let intervals = vec![
            Interval { start: 10, end: 20 },
            Interval { start: 30, end: 40 },
        ];
        assert!(in_targets(&intervals, 35));
    }

    #[test]
    fn in_targets_empty_slice_returns_false() {
        assert!(!in_targets(&[], 0));
        assert!(!in_targets(&[], 999));
    }

    #[test]
    fn load_targets_merges_overlapping_intervals() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "chr1\t10\t20")?;
        writeln!(f, "chr1\t15\t30")?;

        let mut ref_map = HashMap::new();
        ref_map.insert("chr1".to_string(), 0usize);
        let targets = load_targets(f.path(), &ref_map)?;

        let intervals = &targets.by_ref[&0];
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 10);
        assert_eq!(intervals[0].end, 30);
        Ok(())
    }

    #[test]
    fn load_targets_merges_adjacent_intervals() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "chr1\t0\t10")?;
        writeln!(f, "chr1\t10\t20")?; // touching at boundary — merged

        let mut ref_map = HashMap::new();
        ref_map.insert("chr1".to_string(), 0usize);
        let targets = load_targets(f.path(), &ref_map)?;

        let intervals = &targets.by_ref[&0];
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].end, 20);
        Ok(())
    }

    #[test]
    fn load_targets_skips_zero_length_intervals() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "chr1\t10\t10")?; // end == start, invalid
        writeln!(f, "chr1\t20\t30")?;

        let mut ref_map = HashMap::new();
        ref_map.insert("chr1".to_string(), 0usize);
        let targets = load_targets(f.path(), &ref_map)?;

        let intervals = &targets.by_ref[&0];
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 20);
        Ok(())
    }

    #[test]
    fn load_targets_skips_comments_and_blanks() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "# comment")?;
        writeln!(f)?;
        writeln!(f, "chr1\t0\t100")?;

        let mut ref_map = HashMap::new();
        ref_map.insert("chr1".to_string(), 0usize);
        let targets = load_targets(f.path(), &ref_map)?;

        assert_eq!(targets.by_ref[&0].len(), 1);
        Ok(())
    }

    #[test]
    fn load_targets_errors_on_unknown_chrom() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chrX\t0\t100").unwrap();

        let ref_map: HashMap<String, usize> = HashMap::new();
        let err = load_targets(f.path(), &ref_map).unwrap_err();
        assert!(err.to_string().contains("not in BAM header"));
    }
}
