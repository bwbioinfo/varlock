//! Native filter format boundaries. These tests require Cargo to execute.
use std::{fs::File, io::Read, path::Path, process::{Command, Output}};

use anyhow::{Context, Result};
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, BinningIndex};
use noodles_vcf::{self as vcf, variant::io::Write as _};
use tempfile::tempdir;

const INPUT: &str = concat!(
    "##fileformat=VCFv4.3\n",
    "##source=native-filter-fixture\n",
    "##contig=<ID=chr1,length=100000>\n",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n",
    "##INFO=<ID=KEEP,Number=0,Type=Flag,Description=\"Keep\">\n",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal\n",
    "chr1\t10\trs1\tA\tC,G\t12.5\tPASS\tAF=0.125,.;KEEP\tGT:DP:AD\t1|2:12:3,4,5\t0/.:.:.,.,.\n",
    "chr1\t20\t.\tA\tG\t.\tPASS\tAF=0.5\tGT:DP:AD\t0/0:8:8,0\t0/0:10:10,0\n",
    "chr1\t30\t.\tA\t<DEL>\t.\tPASS\tEND=20000;KEEP\tGT:DP:AD\t0/1:15:7,8\t./.:.:.,.\n",
);

fn run(input: &Path, output: &Path, input_format: &str, output_format: &str, extra: &[&str]) -> Result<Output> {
    Ok(Command::new(env!("CARGO_BIN_EXE_varlock"))
        .arg("filter").arg("--input").arg(input).arg("--output").arg(output)
        .args(["--input-format", input_format, "--output-format", output_format])
        .args(extra).output()?)
}

fn success(result: Output) {
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
}

fn failure(result: Output, needle: &str) {
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains(needle), "{}", String::from_utf8_lossy(&result.stderr));
}

fn read_bcf(path: &Path) -> Result<(vcf::Header, String)> {
    let mut reader = bcf::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut writer = vcf::io::Writer::new(Vec::new());
    for result in reader.records() {
        writer.write_variant_record(&header, &result?)?;
    }
    Ok((header, String::from_utf8(writer.into_inner())?))
}

#[test]
fn filter_converts_vcf_bcf_and_preserves_selected_typed_records() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    // Deliberately not .bcf: flags, not extensions, select the codec.
    let encoded = dir.path().join("native.data");
    let selected = dir.path().join("selected.data");
    let decoded = dir.path().join("decoded.vcf.gz");
    std::fs::write(&input, INPUT)?;
    success(run(&input, &encoded, "vcf", "bcf", &[])?);
    let mut magic = [0; 5];
    bgzf::io::Reader::new(File::open(&encoded)?).read_exact(&mut magic)?;
    assert_eq!(&magic, b"BCF\x02\x02");
    let (header, records) = read_bcf(&encoded)?;
    assert_eq!(header.sample_names().iter().map(String::as_str).collect::<Vec<_>>(), ["Tumor", "Normal"]);
    assert!(header.infos().contains_key("KEEP"));
    assert!(header.formats().contains_key("AD"));
    assert!(records.contains("AF=0.125,.;KEEP"));
    assert!(records.contains("1|2:12:3,4,5\t0/.:.:.,.,."));
    success(run(&encoded, &selected, "bcf", "bcf", &["--require-info", "KEEP", "--sample-min-dp", "Tumor=10"])?);
    let (_, selected_records) = read_bcf(&selected)?;
    let expected: String = records.lines().filter(|line| !line.starts_with("chr1\t20\t"))
        .map(|line| format!("{line}\n")).collect();
    assert_eq!(selected_records, expected);
    success(run(&selected, &decoded, "bcf", "vcf", &["--index-type", "tbi"])?);
    let mut text = String::new();
    bgzf::io::Reader::new(File::open(&decoded)?).read_to_string(&mut text)?;
    assert!(text.contains("##source=native-filter-fixture"));
    let decoded_records: String = text.lines().filter(|line| !line.starts_with('#'))
        .map(|line| format!("{line}\n")).collect();
    assert_eq!(decoded_records, expected);
    assert!(dir.path().join("decoded.vcf.gz.tbi").exists());
    // CSI must be a native BCF index, not a Tabix/VCF CSI disguised by suffix.
    let index = csi::fs::read(dir.path().join("selected.data.csi"))?;
    assert!(index.header().is_none());
    let mut reader = bcf::io::Reader::new(File::open(&selected)?);
    let header = reader.read_header()?;
    // Query beyond the POS bin, inside INFO/END, to verify span indexing.
    let region = "chr1:19000-19001".parse()?;
    let mut query = reader.query(&header, &index, &region)?;
    assert_eq!(query.records().collect::<std::io::Result<Vec<_>>>()?.len(), 1);
    Ok(())
}

#[test]
fn filter_native_empty_selection_has_valid_header_and_csi() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    let output = dir.path().join("empty.bcf");
    std::fs::write(&input, INPUT)?;
    success(run(&input, &output, "vcf", "bcf", &["--require-info", "ABSENT"])?);
    let (header, records) = read_bcf(&output)?;
    assert_eq!(header.sample_names().len(), 2);
    assert!(records.is_empty());
    let index = csi::fs::read(dir.path().join("empty.bcf.csi"))?;
    assert!(index.header().is_none());
    assert_eq!(index.reference_sequences().len(), 1);
    Ok(())
}

#[test]
fn filter_native_rejects_tbi_wrong_codec_and_unknown_samples_before_output() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    let output = dir.path().join("output.bcf");
    std::fs::write(&input, INPUT)?;
    failure(run(&input, &output, "vcf", "bcf", &["--index-type", "tbi"])?, "BCF output requires");
    assert!(!output.exists());
    failure(run(&input, &output, "bcf", "bcf", &[])?, "header");
    assert!(!output.exists());
    failure(run(&input, &output, "vcf", "bcf", &["--sample-has-alt", "Absent"])?, "Absent");
    assert!(!output.exists());
    Ok(())
}

#[test]
fn filter_native_refuses_existing_output_and_stale_indexes() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    std::fs::write(&input, INPUT)?;
    failure(run(&input, &input, "vcf", "bcf", &[])?, "must be new paths");
    assert_eq!(std::fs::read_to_string(&input)?, INPUT);
    let output = dir.path().join("output.bcf");
    std::fs::write(dir.path().join("output.bcf.csi"), b"old index")?;
    failure(run(&input, &output, "vcf", "bcf", &[])?, "must be new paths");
    assert!(!output.exists());
    Ok(())
}

#[test]
fn filter_native_rejects_missing_declarations_and_unsorted_output() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    std::fs::write(&input, INPUT.replace("##contig=<ID=chr1,length=100000>\n", ""))?;
    failure(run(&input, &dir.path().join("undeclared.bcf"), "vcf", "bcf", &[])?, "contig declaration");
    let header: String = INPUT.lines().filter(|line| line.starts_with('#')).map(|s| format!("{s}\n")).collect();
    let records: String = INPUT.lines().filter(|line| !line.starts_with('#')).rev().map(|s| format!("{s}\n")).collect();
    std::fs::write(&input, format!("{header}{records}"))?;
    failure(run(&input, &dir.path().join("unsorted.bcf"), "vcf", "bcf", &[])?, "coordinate-sorted");
    assert!(!dir.path().join("unsorted.bcf.csi").exists());
    Ok(())
}

#[test]
fn filter_native_rejects_truncated_bcf_records() -> Result<()> {
    let dir = tempdir()?;
    let input = dir.path().join("input.vcf");
    let encoded = dir.path().join("input.bcf");
    std::fs::write(&input, INPUT)?;
    success(run(&input, &encoded, "vcf", "bcf", &[])?);
    let mut raw = Vec::new();
    bgzf::io::Reader::new(File::open(&encoded)?).read_to_end(&mut raw)?;
    raw.truncate(raw.len().checked_sub(1).context("empty fixture")?);
    let truncated = dir.path().join("truncated.bcf");
    let mut writer = bgzf::io::Writer::new(File::create(&truncated)?);
    std::io::Write::write_all(&mut writer, &raw)?;
    writer.try_finish()?;
    let output = dir.path().join("output.bcf");
    failure(run(&truncated, &output, "bcf", "bcf", &[])?, "invalid input record");
    assert!(!dir.path().join("output.bcf.csi").exists());
    Ok(())
}

#[test]
fn other_commands_do_not_advertise_or_accept_bcf_format_flags() -> Result<()> {
    for command in ["call-targets", "annotate", "intersect"] {
        let output = Command::new(env!("CARGO_BIN_EXE_varlock"))
            .args([command, "--output-format", "bcf"]).output()?;
        failure(output, "--output-format");
    }
    Ok(())
}
