mod support;

use std::{
    fs::{self, File},
    path::Path,
    process::{Command, Output},
};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::alignment::record::cigar::{Op, op::Kind};
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

fn call_targets_command(
    bams: &[&Path],
    reference: &Path,
    targets: &Path,
    output: &Path,
    index_type: &str,
    extra_args: &[&str],
) -> Command {
    call_targets_command_with_backend(
        bams, reference, targets, output, index_type, extra_args, true,
    )
}

fn call_targets_command_with_backend(
    bams: &[&Path],
    reference: &Path,
    targets: &Path,
    output: &Path,
    index_type: &str,
    extra_args: &[&str],
    force_cpu: bool,
) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_varlock"));
    command.arg("call-targets");
    #[cfg(feature = "wgpu")]
    if force_cpu {
        command.arg("--cpu");
    }
    #[cfg(not(feature = "wgpu"))]
    let _ = force_cpu;
    for bam in bams {
        command.arg("--input").arg(bam);
    }
    command
        .arg("--reference")
        .arg(reference)
        .arg("--targets")
        .arg(targets)
        .arg("--output")
        .arg(output)
        .arg("--index-type")
        .arg(index_type);
    command.args(extra_args);
    command
}

fn run_call_targets(
    bams: &[&Path],
    reference: &Path,
    targets: &Path,
    output: &Path,
    index_type: &str,
    extra_args: &[&str],
) -> Result<Output> {
    let result = call_targets_command(bams, reference, targets, output, index_type, extra_args)
        .output()
        .context("failed to run varlock call-targets")?;

    assert!(result.status.success(), "{}", output_text(&result));
    Ok(result)
}

fn sample_names(parsed: &support::call_targets::ParsedVcf) -> Vec<&str> {
    parsed
        .header
        .sample_names()
        .iter()
        .map(String::as_str)
        .collect()
}

fn sample_depth(parsed: &support::call_targets::ParsedVcf, sample_name: &str) -> Result<i32> {
    let [record] = parsed.records.as_slice() else {
        anyhow::bail!("expected one VCF record, got {}", parsed.records.len());
    };
    let sample = record
        .samples()
        .get(&parsed.header, sample_name)
        .with_context(|| format!("missing {sample_name} values"))?;
    match sample.get("DP") {
        Some(Some(SampleValue::Integer(depth))) => Ok(*depth),
        value => anyhow::bail!("expected integer DP for {sample_name}, got {value:?}"),
    }
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
        run_call_targets(
            &[&bam.path],
            &reference.path,
            &targets,
            &output,
            index_type,
            &[],
        )?;

        let parsed = read_vcf(&output)?;
        assert_eq!(sample_names(&parsed), ["sample_a"]);

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
fn call_targets_emits_indels_by_default_and_honors_indel_mode_flags() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"ACGTACGTACGTACGTACGT")],
    )?;
    // The fifth base is the VCF anchor for all three synthetic calls.
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("rg0", Some("sample_a"));
    builder.add_read(ReadSpec::single_base("snv", 5, b'C').with_read_group("rg0"));
    builder.add_read(
        ReadSpec::single_base("insertion", 5, b'A')
            .with_sequence(b"AGC")
            .with_cigar(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Insertion, 1),
                Op::new(Kind::Match, 1),
            ])
            .with_read_group("rg0"),
    );
    builder.add_read(
        ReadSpec::single_base("deletion", 5, b'A')
            .with_sequence(b"AG")
            .with_cigar(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 1),
            ])
            .with_read_group("rg0"),
    );
    let bam = builder.write(dir.path().join("sample.bam"))?;

    let default_output = dir.path().join("default.vcf.gz");
    run_call_targets(
        &[&bam.path],
        &reference.path,
        &targets,
        &default_output,
        "csi",
        &[],
    )?;
    let default_calls = read_vcf(&default_output)?;
    let alleles = default_calls
        .records
        .iter()
        .map(|record| {
            (
                record.variant_start().unwrap().get(),
                record.reference_bases().to_string(),
                record.alternate_bases().as_ref()[0].clone(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        alleles,
        [
            (5, "A".to_string(), "C".to_string()),
            (5, "A".to_string(), "AG".to_string()),
            (5, "AC".to_string(), "A".to_string())
        ]
    );
    for record in &default_calls.records {
        let sample = record
            .samples()
            .get(&default_calls.header, "sample_a")
            .context("missing sample_a values")?;
        assert_eq!(sample.get("DP"), Some(Some(&SampleValue::Integer(3))));
        assert_eq!(
            sample.get("AD"),
            Some(Some(&SampleValue::Array(SampleArray::Integer(vec![
                Some(2),
                Some(1),
            ]))))
        );
    }

    #[cfg(feature = "wgpu")]
    {
        let gpu_output = dir.path().join("gpu.vcf.gz");
        let result = call_targets_command_with_backend(
            &[&bam.path],
            &reference.path,
            &targets,
            &gpu_output,
            "csi",
            &[],
            false,
        )
        .output()?;
        assert!(result.status.success(), "{}", output_text(&result));
        let gpu_calls = read_vcf(&gpu_output)?;
        let gpu_alleles = gpu_calls
            .records
            .iter()
            .map(|record| {
                (
                    record.variant_start().unwrap().get(),
                    record.reference_bases().to_string(),
                    record.alternate_bases().as_ref()[0].clone(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(gpu_alleles, alleles);
    }

    let no_indels_output = dir.path().join("no-indels.vcf.gz");
    run_call_targets(
        &[&bam.path],
        &reference.path,
        &targets,
        &no_indels_output,
        "csi",
        &["--no-indels"],
    )?;
    let no_indels = read_vcf(&no_indels_output)?;
    assert_eq!(no_indels.records.len(), 1);
    assert_eq!(no_indels.records[0].alternate_bases().as_ref()[0], "C");

    let indels_only_output = dir.path().join("indels-only.vcf.gz");
    run_call_targets(
        &[&bam.path],
        &reference.path,
        &targets,
        &indels_only_output,
        "csi",
        &["--indels-only"],
    )?;
    let indels_only = read_vcf(&indels_only_output)?;
    assert_eq!(indels_only.records.len(), 2);
    assert!(
        indels_only
            .records
            .iter()
            .all(|record| record.alternate_bases().as_ref()[0] != "C")
    );

    let conflicting_output = dir.path().join("conflicting.vcf.gz");
    let result = call_targets_command(
        &[&bam.path],
        &reference.path,
        &targets,
        &conflicting_output,
        "csi",
        &["--no-indels", "--indels-only"],
    )
    .output()?;
    assert!(!result.status.success(), "{}", output_text(&result));
    assert!(output_text(&result).contains("cannot be used with"));

    Ok(())
}

#[test]
fn call_targets_split_by_preserves_file_and_read_group_identity_without_sm() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut first = BamFixtureBuilder::new("chr1", 20);
    first.add_read_group("rg1", None);
    first.add_read_group("rg2", None);
    first.add_observations(5, b'C', 1, Some("rg1"));
    first.add_observations(5, b'C', 1, Some("rg2"));
    let first = first.write(dir.path().join("first.bam"))?;

    let mut second = BamFixtureBuilder::new("chr1", 20);
    second.add_read_group("rg1", None);
    second.add_observations(5, b'C', 1, Some("rg1"));
    let second = second.write(dir.path().join("second.bam"))?;

    let policies: &[(&str, &[&str], &[i32])] = &[
        ("sm", &["SAMPLE"], &[3]),
        ("file", &["first", "second"], &[2, 1]),
        (
            "rg",
            &["first__rg1", "first__rg2", "second__rg1"],
            &[1, 1, 1],
        ),
    ];
    for index_type in ["csi", "tbi"] {
        for &(split_by, expected_names, expected_depths) in policies {
            let output = dir
                .path()
                .join(format!("calls.{split_by}.{index_type}.vcf.gz"));
            run_call_targets(
                &[&first.path, &second.path],
                &reference.path,
                &targets,
                &output,
                index_type,
                &["--split-by", split_by],
            )?;
            let parsed = read_vcf(&output)?;
            assert_eq!(sample_names(&parsed), expected_names);
            for (&sample_name, &expected_depth) in expected_names.iter().zip(expected_depths) {
                assert_eq!(sample_depth(&parsed, sample_name)?, expected_depth);
            }
            assert!(append_extension(&output, index_type).exists());
        }
    }

    Ok(())
}

#[test]
fn call_targets_sm_scopes_reused_read_groups_per_input_and_merges_equal_sm() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut first = BamFixtureBuilder::new("chr1", 20);
    first.add_read_group("shared", Some("sample_a"));
    first.add_observations(5, b'C', 1, Some("shared"));
    let first = first.write(dir.path().join("first.bam"))?;

    let mut second = BamFixtureBuilder::new("chr1", 20);
    second.add_read_group("shared", Some("sample_b"));
    second.add_observations(5, b'C', 1, Some("shared"));
    let second = second.write(dir.path().join("second.bam"))?;

    let mut third = BamFixtureBuilder::new("chr1", 20);
    third.add_read_group("other", Some("sample_a"));
    third.add_observations(5, b'C', 1, Some("other"));
    let third = third.write(dir.path().join("third.bam"))?;

    let output = dir.path().join("calls.vcf.gz");
    run_call_targets(
        &[&first.path, &second.path, &third.path],
        &reference.path,
        &targets,
        &output,
        "csi",
        &[],
    )?;
    let parsed = read_vcf(&output)?;
    assert_eq!(sample_names(&parsed), ["sample_a", "sample_b"]);
    assert_eq!(sample_depth(&parsed, "sample_a")?, 2);
    assert_eq!(sample_depth(&parsed, "sample_b")?, 1);

    Ok(())
}

#[test]
fn call_targets_sm_keeps_inputs_without_sm_when_other_inputs_have_sm() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut mapped = BamFixtureBuilder::new("chr1", 20);
    mapped.add_read_group("mapped", Some("sample_a"));
    mapped.add_observations(5, b'C', 1, Some("mapped"));
    let mapped = mapped.write(dir.path().join("mapped.bam"))?;

    let mut unmapped = BamFixtureBuilder::new("chr1", 20);
    unmapped.add_read_group("unmapped", None);
    unmapped.add_observations(5, b'C', 1, Some("unmapped"));
    let unmapped = unmapped.write(dir.path().join("unmapped.bam"))?;

    let output = dir.path().join("calls.vcf.gz");
    run_call_targets(
        &[&mapped.path, &unmapped.path],
        &reference.path,
        &targets,
        &output,
        "csi",
        &[],
    )?;
    let parsed = read_vcf(&output)?;
    assert_eq!(sample_names(&parsed), ["sample_a", "unmapped"]);
    assert_eq!(sample_depth(&parsed, "sample_a")?, 1);
    assert_eq!(sample_depth(&parsed, "unmapped")?, 1);

    Ok(())
}

#[test]
fn call_targets_accepts_rg_map_only_with_sm_split_by() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;
    let mut builder = BamFixtureBuilder::new("chr1", 20);
    builder.add_read_group("rg1", None);
    builder.add_observations(5, b'C', 1, Some("rg1"));
    let bam = builder.write(dir.path().join("sample.bam"))?;
    let map = dir.path().join("rg-map.tsv");
    fs::write(&map, "RG\tSM\nrg1\tsample_a\n")?;

    for split_by in ["file", "rg"] {
        let output = dir.path().join(format!("calls.{split_by}.vcf.gz"));
        let result = call_targets_command(
            &[&bam.path],
            &reference.path,
            &targets,
            &output,
            "csi",
            &["--split-by", split_by, "--rg-map", map.to_str().unwrap()],
        )
        .output()?;
        assert!(!result.status.success(), "{}", output_text(&result));
        assert!(output_text(&result).contains("--rg-map is only supported with --split-by sm"));
    }

    let output = dir.path().join("calls.sm.vcf.gz");
    run_call_targets(
        &[&bam.path],
        &reference.path,
        &targets,
        &output,
        "csi",
        &[
            "--split-by",
            "sm",
            "--rg-map",
            map.to_str().context("invalid map path")?,
        ],
    )?;
    let parsed = read_vcf(&output)?;
    assert_eq!(sample_names(&parsed), ["sample_a"]);
    assert_eq!(sample_depth(&parsed, "sample_a")?, 1);

    Ok(())
}

#[test]
fn call_targets_rg_split_rejects_missing_or_unknown_read_groups() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut missing = BamFixtureBuilder::new("chr1", 20);
    missing.add_read_group("declared", None);
    missing.add_observations(5, b'C', 1, None);
    let missing = missing.write(dir.path().join("missing.bam"))?;

    let mut unknown = BamFixtureBuilder::new("chr1", 20);
    unknown.add_read_group("declared", None);
    unknown.add_observations(5, b'C', 1, Some("unknown"));
    let unknown = unknown.write(dir.path().join("unknown.bam"))?;

    for (bam, expected_error) in [
        (&missing, "requires an RG tag"),
        (&unknown, "is not declared in its BAM header"),
    ] {
        let output = dir.path().join(format!(
            "{}.vcf.gz",
            bam.path.file_stem().unwrap().display()
        ));
        let result = call_targets_command(
            &[&bam.path],
            &reference.path,
            &targets,
            &output,
            "csi",
            &["--split-by", "rg"],
        )
        .output()?;
        assert!(!result.status.success(), "{}", output_text(&result));
        assert!(output_text(&result).contains(expected_error));
        assert!(!output.exists());
    }

    Ok(())
}

#[cfg(feature = "wgpu")]
#[test]
fn call_targets_gpu_static_path_matches_cpu_sample_resolution() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_fasta(
        dir.path().join("reference.fa"),
        &[("chr1", b"AAAAAAAAAAAAAAAAAAAA")],
    )?;
    let targets = write_bed(dir.path().join("targets.bed"), &[("chr1", 4, 5)])?;

    let mut first = BamFixtureBuilder::new("chr1", 20);
    first.add_read_group("shared", None);
    first.add_observations(5, b'C', 1, Some("shared"));
    let first = first.write(dir.path().join("first.bam"))?;

    let mut second = BamFixtureBuilder::new("chr1", 20);
    second.add_read_group("shared", None);
    second.add_observations(5, b'C', 1, Some("shared"));
    let second = second.write(dir.path().join("second.bam"))?;

    let cpu_output = dir.path().join("cpu.vcf.gz");
    run_call_targets(
        &[&first.path, &second.path],
        &reference.path,
        &targets,
        &cpu_output,
        "csi",
        &["--split-by", "rg"],
    )?;

    let gpu_output = dir.path().join("gpu.vcf.gz");
    let result = call_targets_command_with_backend(
        &[&first.path, &second.path],
        &reference.path,
        &targets,
        &gpu_output,
        "csi",
        &["--split-by", "rg"],
        false,
    )
    .output()?;
    assert!(result.status.success(), "{}", output_text(&result));

    let cpu = read_vcf(&cpu_output)?;
    let gpu = read_vcf(&gpu_output)?;
    assert_eq!(sample_names(&gpu), sample_names(&cpu));
    for sample_name in sample_names(&cpu) {
        assert_eq!(
            sample_depth(&gpu, sample_name)?,
            sample_depth(&cpu, sample_name)?
        );
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
