#[allow(dead_code)]
mod support;

use std::{
    fs::File,
    path::Path,
    process::{Command, Output},
};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::header::record::value::map::read_group::tag::SAMPLE;
use tempfile::tempdir;

use support::call_targets::BamFixtureBuilder;

fn run_varlock(args: &[&str]) -> Result<Output> {
    Command::new(env!("CARGO_BIN_EXE_varlock"))
        .args(args)
        .output()
        .context("failed to run varlock")
}

fn output_text(output: &Output) -> String {
    format!(
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn sample_tag<'a>(header: &'a noodles_sam::Header, read_group: &str) -> Option<&'a [u8]> {
    header
        .read_groups()
        .get(read_group.as_bytes())
        .and_then(|group| group.other_fields().get(&SAMPLE))
        .map(|value| value.as_ref())
}

fn record_count(path: &Path) -> Result<usize> {
    let mut reader = bam::io::Reader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    reader.read_header()?;
    reader
        .records()
        .collect::<std::io::Result<Vec<_>>>()
        .map(|records| records.len())
        .map_err(Into::into)
}

fn observed_bases(path: &Path) -> Result<Vec<u8>> {
    let mut reader = bam::io::Reader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    reader.read_header()?;
    reader
        .records()
        .map(|result| result.map(|record| record.sequence().get(0).unwrap_or(b'N')))
        .collect::<std::io::Result<Vec<_>>>()
        .map_err(Into::into)
}

#[test]
fn add_sm_to_bam_defaults_to_input_stem_and_reindexes_output() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("cohort-01.bam");
    let output = dir.path().join("rewritten.bam");

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("missing-a", None);
    builder.add_read_group("missing-b", None);
    builder.add_read_group("existing", Some("already-set"));
    builder.add_observations(5, b'A', 1, Some("missing-a"));
    builder.add_observations(5, b'C', 1, Some("missing-b"));
    builder.add_observations(5, b'G', 1, Some("existing"));
    let fixture = builder.write(&input)?;

    let result = run_varlock(&[
        "add-sm-to-bam",
        "--input",
        fixture.path.to_str().context("invalid input path")?,
        "--output",
        output.to_str().context("invalid output path")?,
    ])?;
    assert!(result.status.success(), "{}", output_text(&result));
    assert!(output.exists());
    assert!(output.with_extension("bam.bai").exists());

    let mut reader = bam::io::Reader::new(File::open(&output)?);
    let header = reader.read_header()?;
    assert_eq!(sample_tag(&header, "missing-a"), Some(&b"cohort-01"[..]));
    assert_eq!(sample_tag(&header, "missing-b"), Some(&b"cohort-01"[..]));
    assert_eq!(sample_tag(&header, "existing"), Some(&b"already-set"[..]));
    assert_eq!(record_count(&output)?, 3);
    assert_eq!(observed_bases(&output)?, observed_bases(&fixture.path)?);

    let mut indexed_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&output)
        .context("output BAI should be readable")?;
    assert_eq!(indexed_reader.read_header()?.read_groups().len(), 3);

    Ok(())
}

#[test]
fn add_sm_to_bam_refuses_to_overwrite_the_input() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("sample.bam");
    let output_alias = dir.path().join(".").join("sample.bam");

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("rg0", None);
    builder.add_observations(5, b'A', 1, Some("rg0"));
    let fixture = builder.write(&input)?;

    let result = run_varlock(&[
        "add-sm-to-bam",
        "--input",
        fixture.path.to_str().context("invalid input path")?,
        "--output",
        output_alias.to_str().context("invalid output alias")?,
    ])?;
    assert!(!result.status.success(), "{}", output_text(&result));
    assert!(output_text(&result).contains("--output must differ from --input"));

    let mut reader = bam::io::Reader::new(File::open(&fixture.path)?);
    let header = reader.read_header()?;
    assert_eq!(sample_tag(&header, "rg0"), None);

    Ok(())
}

#[test]
fn add_sm_to_bam_uses_explicit_sample_for_missing_tags() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("source.bam");
    let output = dir.path().join("output.bam");

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("rg0", None);
    builder.add_observations(5, b'A', 1, Some("rg0"));
    let fixture = builder.write(&input)?;

    let result = run_varlock(&[
        "add-sm-to-bam",
        "--input",
        fixture.path.to_str().context("invalid input path")?,
        "--output",
        output.to_str().context("invalid output path")?,
        "--sample",
        "override",
    ])?;
    assert!(result.status.success(), "{}", output_text(&result));

    let mut reader = bam::io::Reader::new(File::open(&output)?);
    let header = reader.read_header()?;
    assert_eq!(sample_tag(&header, "rg0"), Some(&b"override"[..]));

    Ok(())
}

#[test]
fn add_sm_to_bam_rejects_bams_without_read_groups_before_writing_output() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("no-read-groups.bam");
    let output = dir.path().join("output.bam");

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_observations(5, b'A', 1, None);
    let fixture = builder.write(&input)?;

    let result = run_varlock(&[
        "add-sm-to-bam",
        "--input",
        fixture.path.to_str().context("invalid input path")?,
        "--output",
        output.to_str().context("invalid output path")?,
    ])?;
    assert!(!result.status.success(), "{}", output_text(&result));
    assert!(output_text(&result).contains("does not contain any @RG records"));
    assert!(!output.exists());

    Ok(())
}
