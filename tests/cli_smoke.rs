use std::{
    fs::File,
    io::Read,
    process::{Command, Output},
};

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use tempfile::tempdir;

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

#[test]
fn top_level_help_lists_core_commands() -> Result<()> {
    let output = run_varlock(&["--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("call-targets"));
    assert!(stdout.contains("annotate"));
    assert!(stdout.contains("filter"));
    assert!(stdout.contains("intersect"));
    Ok(())
}

#[test]
fn annotate_help_lists_index_type() -> Result<()> {
    let output = run_varlock(&["annotate", "--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--reference"));
    assert!(stdout.contains("--index-type"));
    assert!(stdout.contains("possible values: csi, tbi"));
    Ok(())
}

#[test]
fn call_targets_help_lists_pair_flags() -> Result<()> {
    let output = run_varlock(&["call-targets", "--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--pair"));
    assert!(stdout.contains("--tumor-min-alt-count"));
    assert!(stdout.contains("--normal-max-alt-fraction"));
    Ok(())
}

#[cfg(feature = "wgpu")]
#[test]
fn call_targets_gpu_help_lists_pair_flags() -> Result<()> {
    let output = run_varlock(&["call-targets-gpu", "--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--pair"));
    assert!(stdout.contains("--tumor-min-alt-count"));
    assert!(stdout.contains("--normal-max-alt-fraction"));
    Ok(())
}

#[test]
fn filter_help_lists_info_predicates() -> Result<()> {
    let output = run_varlock(&["filter", "--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--require-info"));
    assert!(stdout.contains("--exclude-info"));
    assert!(stdout.contains("--max-info"));
    assert!(stdout.contains("--expr"));
    assert!(stdout.contains("--sample-group"));
    assert!(stdout.contains("--sample-has-alt"));
    assert!(stdout.contains("--group-all-min-dp"));
    Ok(())
}

#[test]
fn intersect_help_lists_two_file_modes() -> Result<()> {
    let output = run_varlock(&["intersect", "--help"])?;
    assert!(output.status.success(), "{}", output_text(&output));

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--left"));
    assert!(stdout.contains("--right"));
    assert!(stdout.contains("--input"));
    assert!(stdout.contains("--left-samples"));
    assert!(stdout.contains("--right-samples"));
    assert!(stdout.contains("--mode"));
    assert!(stdout.contains("possible values: shared, left-only, right-only"));
    Ok(())
}

#[test]
fn annotate_smoke_writes_bgzipped_vcf_and_index() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    let database = dir.path().join("db.vcf");
    let output = dir.path().join("out.vcf.gz");

    std::fs::write(
        &input,
        "##fileformat=VCFv4.3\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t10\t.\tA\tC\t.\tPASS\t.\n",
    )?;
    std::fs::write(
        &database,
        "##fileformat=VCFv4.3\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25\n",
    )?;

    let output_result = run_varlock(&[
        "annotate",
        "--input",
        input.to_str().context("invalid input path")?,
        "--database",
        &format!("db={}", database.display()),
        "--annotation",
        "db:AF=db_AF",
        "--output",
        output.to_str().context("invalid output path")?,
    ])?;
    assert!(
        output_result.status.success(),
        "{}",
        output_text(&output_result)
    );

    assert!(output.exists());
    assert!(dir.path().join("out.vcf.gz.csi").exists());

    let mut reader = bgzf::io::Reader::new(File::open(output)?);
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    assert!(text.contains("##INFO=<ID=db_AF,Number=A"));
    assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tdb_AF=0.25"));
    Ok(())
}

#[test]
fn filter_smoke_writes_filtered_bgzipped_vcf_and_index() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    let output = dir.path().join("out.vcf.gz");

    std::fs::write(
        &input,
        "##fileformat=VCFv4.3\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05\n\
         chr1\t11\t.\tA\tG\t.\tPASS\tAF=0.20\n",
    )?;

    let output_result = run_varlock(&[
        "filter",
        "--input",
        input.to_str().context("invalid input path")?,
        "--max-info",
        "AF=0.10",
        "--output",
        output.to_str().context("invalid output path")?,
    ])?;
    assert!(
        output_result.status.success(),
        "{}",
        output_text(&output_result)
    );

    assert!(output.exists());
    assert!(dir.path().join("out.vcf.gz.csi").exists());

    let mut reader = bgzf::io::Reader::new(File::open(output)?);
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05"));
    assert!(!text.contains("chr1\t11\t.\tA\tG"));
    Ok(())
}
