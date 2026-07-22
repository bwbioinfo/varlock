mod support;

use std::{fs::File, process::Command};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::header::record::value::map::read_group::tag::SAMPLE;
use noodles_vcf::variant::record_buf::{
    info::field::Value as InfoValue,
    samples::sample::{Value as SampleValue, value::Array as SampleArray},
};
use tempfile::tempdir;

use support::call_targets::{
    BamFixtureBuilder, ReadSpec, append_extension, read_vcf, write_bed, write_fasta,
};

fn output_text(output: &std::process::Output) -> String {
    format!(
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn run_call_targets(
    bam: &std::path::Path,
    reference: &std::path::Path,
    targets: &std::path::Path,
    output: &std::path::Path,
    index_type: &str,
) -> Result<()> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_varlock"));
    command.arg("call-targets");
    #[cfg(feature = "wgpu")]
    command.arg("--cpu");
    let result = command
        .arg("--input")
        .arg(bam)
        .arg("--reference")
        .arg(reference)
        .arg("--targets")
        .arg(targets)
        .arg("--output")
        .arg(output)
        .arg("--index-type")
        .arg(index_type)
        .output()
        .context("failed to run varlock call-targets")?;

    assert!(result.status.success(), "{}", output_text(&result));
    Ok(())
}

#[test]
fn call_targets_cpu_reads_bam_and_writes_typed_indexed_vcf() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("rg0", Some("sample_a"));
    builder.add_observations(5, b'A', 2, Some("rg0"));
    builder.add_observations(5, b'C', 3, Some("rg0"));
    let bam = builder.write(dir.path().join("sample.bam"))?;

    assert!(reference.fai_path.exists());
    assert!(bam.bai_path.exists());

    let mut indexed_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&bam.path)
        .context("failed to open generated BAM and BAI")?;
    let bam_header = indexed_reader.read_header()?;
    let region: Region = "chr1:5-5".parse()?;
    let query = indexed_reader.query(&bam_header, &region)?;
    assert_eq!(
        query.records().collect::<std::io::Result<Vec<_>>>()?.len(),
        5
    );

    for index_type in ["csi", "tbi"] {
        let output = dir.path().join(format!("calls.{index_type}.vcf.gz"));
        run_call_targets(&bam.path, &reference.path, &targets, &output, index_type)?;

        let parsed = read_vcf(&output)?;
        let sample_names = parsed
            .header
            .sample_names()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        assert_eq!(sample_names, ["sample_a"]);

        let [record] = parsed.records.as_slice() else {
            panic!("expected one VCF record, got {}", parsed.records.len());
        };
        assert_eq!(record.reference_sequence_name(), "chr1");
        assert_eq!(record.variant_start(), Some(Position::try_from(5)?));
        assert_eq!(record.reference_bases(), "A");
        assert_eq!(
            record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["C"]
        );
        assert_eq!(record.info().get("DP"), Some(Some(&InfoValue::Integer(5))));

        let sample = record
            .samples()
            .get(&parsed.header, "sample_a")
            .context("missing sample_a values")?;
        assert_eq!(
            sample.get("GT"),
            Some(Some(&SampleValue::Genotype("0/1".parse()?)))
        );
        assert_eq!(sample.get("DP"), Some(Some(&SampleValue::Integer(5))));
        assert_eq!(
            sample.get("AD"),
            Some(Some(&SampleValue::Array(SampleArray::Integer(vec![
                Some(2),
                Some(3),
            ]))))
        );

        let index_path = append_extension(&output, index_type);
        match index_type {
            "csi" => {
                let index = noodles_csi::fs::read(&index_path)?;
                assert_eq!(index.reference_sequences().len(), 1);
            }
            "tbi" => {
                let index = noodles_tabix::fs::read(&index_path)?;
                assert_eq!(index.reference_sequences().len(), 1);
            }
            _ => unreachable!(),
        }
    }

    Ok(())
}

#[test]
fn fixture_builder_generates_follow_up_edge_case_shapes() -> Result<()> {
    let dir = tempdir()?;

    let mut first = BamFixtureBuilder::new("chr1", 20);
    first.add_read_group("shared-rg", Some("sample_a"));
    for (i, base) in (*b"CACACA").into_iter().enumerate() {
        first.add_read(
            ReadSpec::single_base(format!("ordered-{i}"), 5, base).with_read_group("shared-rg"),
        );
    }
    let first = first.write(dir.path().join("first.bam"))?;

    let mut duplicate_rg = BamFixtureBuilder::new("chr1", 20);
    duplicate_rg.add_read_group("shared-rg", Some("sample_b"));
    duplicate_rg.add_observations(5, b'C', 1, Some("shared-rg"));
    let duplicate_rg = duplicate_rg.write(dir.path().join("duplicate-rg.bam"))?;

    let mut rgless = BamFixtureBuilder::new("chr1", 20);
    rgless.add_observations(5, b'A', 1, None);
    let rgless = rgless.write(dir.path().join("rgless.bam"))?;

    let mut mismatched_dictionary = BamFixtureBuilder::new("chr1", 19);
    mismatched_dictionary.add_observations(5, b'A', 1, None);
    let mismatched_dictionary =
        mismatched_dictionary.write(dir.path().join("mismatched-dictionary.bam"))?;

    let mut first_reader = bam::io::Reader::new(File::open(&first.path)?);
    let first_header = first_reader.read_header()?;
    let first_rg = first_header
        .read_groups()
        .get(&b"shared-rg"[..])
        .context("missing shared-rg in first BAM")?;
    assert_eq!(
        first_rg
            .other_fields()
            .get(&SAMPLE)
            .map(|value| value.as_ref()),
        Some(&b"sample_a"[..])
    );
    let observed_order = first_reader
        .records()
        .map(|result| result.map(|record| record.sequence().get(0).unwrap_or(b'N')))
        .collect::<std::io::Result<Vec<_>>>()?;
    assert_eq!(observed_order, b"CACACA");
    assert!(
        observed_order.len() > 5,
        "fixture must exceed a depth cap of 5"
    );

    let mut duplicate_reader = bam::io::Reader::new(File::open(&duplicate_rg.path)?);
    let duplicate_header = duplicate_reader.read_header()?;
    let duplicate = duplicate_header
        .read_groups()
        .get(&b"shared-rg"[..])
        .context("missing duplicate shared-rg")?;
    assert_eq!(
        duplicate
            .other_fields()
            .get(&SAMPLE)
            .map(|value| value.as_ref()),
        Some(&b"sample_b"[..])
    );

    let mut rgless_reader = bam::io::Reader::new(File::open(&rgless.path)?);
    assert!(rgless_reader.read_header()?.read_groups().is_empty());

    let mut mismatch_reader = bam::io::Reader::new(File::open(&mismatched_dictionary.path)?);
    let mismatch_header = mismatch_reader.read_header()?;
    let chr1 = mismatch_header
        .reference_sequences()
        .get(&b"chr1"[..])
        .context("missing chr1 from mismatched dictionary")?;
    assert_eq!(usize::from(chr1.length()), 19);

    for fixture in [first, duplicate_rg, rgless, mismatched_dictionary] {
        bam::io::indexed_reader::Builder::default()
            .build_from_path(&fixture.path)
            .with_context(|| format!("fixture index is unreadable: {}", fixture.path.display()))?;
    }

    Ok(())
}
