use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use noodles_bam::io::Reader;
use noodles_sam as sam;

pub(crate) fn collect_samples(
    paths: &[PathBuf],
    rg_map: Option<&[(String, String)]>,
) -> Result<(Vec<String>, HashMap<String, String>)> {
    if let Some(rows) = rg_map {
        return build_rg_to_sm_from_map(rows);
    }

    let mut rg_to_sm = HashMap::new();
    for path in paths {
        let file = File::open(path)
            .with_context(|| format!("failed to open input BAM {}", path.display()))?;
        let mut reader = Reader::new(file);
        let input_header = reader
            .read_header()
            .with_context(|| format!("failed to read header for {}", path.display()))?;
        let from_header = rg_to_sm_from_header(&input_header)?;
        for (rg, sm) in from_header {
            insert_rg_sample(
                &mut rg_to_sm,
                rg,
                sm,
                &format!("BAM header {}", path.display()),
            )?;
        }
    }

    let mut samples: Vec<String> = rg_to_sm.values().cloned().collect();
    samples.sort();
    samples.dedup();

    Ok((samples, rg_to_sm))
}

pub(crate) fn read_rg_map(path: &Path) -> Result<Vec<(String, String)>> {
    let file =
        File::open(path).with_context(|| format!("failed to open rg_map {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line.context("failed to read rg_map line")?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        rows.push((parts[0].to_string(), parts[1].to_string()));
    }
    Ok(rows)
}

fn build_rg_to_sm_from_map(
    rows: &[(String, String)],
) -> Result<(Vec<String>, HashMap<String, String>)> {
    let mut rg_to_sm = HashMap::new();
    let mut samples = Vec::new();

    for (rg, sm) in rows {
        insert_rg_sample(&mut rg_to_sm, rg.clone(), sm.clone(), "rg-map")?;
        if !samples.contains(sm) {
            samples.push(sm.clone());
        }
    }

    if rg_to_sm.is_empty() {
        bail!("rg-map did not contain any RG-to-sample mappings");
    }

    Ok((samples, rg_to_sm))
}

fn insert_rg_sample(
    rg_to_sm: &mut HashMap<String, String>,
    rg: String,
    sm: String,
    source: &str,
) -> Result<()> {
    if let Some(existing) = rg_to_sm.get(&rg) {
        if existing != &sm {
            bail!(
                "conflicting sample mapping for RG {}: {} vs {} ({})",
                rg,
                existing,
                sm,
                source
            );
        }
        return Ok(());
    }
    rg_to_sm.insert(rg, sm);
    Ok(())
}

fn rg_to_sm_from_header(header: &sam::Header) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for (rg_id, read_group) in header.read_groups() {
        let rg_id = std::str::from_utf8(rg_id.as_ref()).context("invalid RG id")?;
        let sm = read_group
            .other_fields()
            .get(&sam::header::record::value::map::read_group::tag::SAMPLE)
            .map(|value| std::str::from_utf8(value.as_ref()).context("invalid SM value"))
            .transpose()?;
        if let Some(sm) = sm {
            map.insert(rg_id.to_string(), sm.to_string());
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use anyhow::Result;
    use tempfile::NamedTempFile;
    use super::{build_rg_to_sm_from_map, read_rg_map};

    #[test]
    fn read_rg_map_parses_tab_separated_lines() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "rg2\tsample_b")?;

        let rows = read_rg_map(f.path())?;
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0], ("rg1".to_string(), "sample_a".to_string()));
        assert_eq!(rows[1], ("rg2".to_string(), "sample_b".to_string()));
        Ok(())
    }

    #[test]
    fn read_rg_map_skips_comments_and_blank_lines() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "# header comment")?;
        writeln!(f)?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "rg2\tsample_b")?;

        let rows = read_rg_map(f.path())?;
        assert_eq!(rows.len(), 2);
        Ok(())
    }

    #[test]
    fn read_rg_map_skips_lines_with_fewer_than_two_fields() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "orphan_rg")?; // only one field — skip
        writeln!(f, "rg2\tsample_b")?;

        let rows = read_rg_map(f.path())?;
        assert_eq!(rows.len(), 2);
        Ok(())
    }

    #[test]
    fn build_rg_to_sm_preserves_insertion_order_for_samples() -> Result<()> {
        let rows = vec![
            ("rg1".to_string(), "sample_a".to_string()),
            ("rg2".to_string(), "sample_b".to_string()),
            ("rg3".to_string(), "sample_b".to_string()), // second RG for same sample
        ];
        let (samples, rg_to_sm) = build_rg_to_sm_from_map(&rows)?;
        assert_eq!(samples, vec!["sample_a", "sample_b"]);
        assert_eq!(rg_to_sm["rg1"], "sample_a");
        assert_eq!(rg_to_sm["rg3"], "sample_b");
        Ok(())
    }

    #[test]
    fn build_rg_to_sm_errors_on_conflicting_mapping() {
        let rows = vec![
            ("rg1".to_string(), "sample_a".to_string()),
            ("rg1".to_string(), "sample_b".to_string()), // same RG, different sample
        ];
        let err = build_rg_to_sm_from_map(&rows).unwrap_err();
        assert!(err.to_string().contains("conflicting sample mapping"));
    }

    #[test]
    fn build_rg_to_sm_errors_on_empty_input() {
        let rows: Vec<(String, String)> = vec![];
        let err = build_rg_to_sm_from_map(&rows).unwrap_err();
        assert!(err.to_string().contains("did not contain any"));
    }
}
