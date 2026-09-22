use std::{
    fs::{self, File},
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, Output},
};

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use tempfile::tempdir;

const HEADER: &str = "##fileformat=VCFv4.3\n##source=contract-fixture\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
const RECORDS: &str = "chr1\t10\trs10\tA\tC\t.\tPASS\tAF=0.05\nchr1\t11\trs11\tA\tG\t.\tPASS\tAF=0.20\n";

fn run(args: &[&str]) -> Result<Output> {
    Command::new(env!("CARGO_BIN_EXE_varlock"))
        .args(args)
        .output()
        .context("failed to run varlock")
}

fn write_plain(path: &Path, body: &str) -> Result<()> {
    fs::write(path, body).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

fn write_bgzf(path: &Path, body: &str) -> Result<()> {
    let mut writer = bgzf::io::Writer::new(File::create(path)?);
    writer.write_all(body.as_bytes())?;
    writer.try_finish()?;
    Ok(())
}

fn read_bgzf(path: &Path) -> Result<String> {
    let mut reader = bgzf::io::Reader::new(File::open(path)?);
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    Ok(text)
}

fn write_reference(dir: &Path) -> Result<PathBuf> {
    let path = dir.join("reference.fa");
    write_plain(&path, ">chr1\nACGTACGTACGTACGTACGT\n")?;
    // FASTA index fields: name, length, sequence byte offset, bases/line, bytes/line.
    write_plain(&PathBuf::from(format!("{}.fai", path.display())), "chr1\t20\t6\t20\t21\n")?;
    Ok(path)
}

fn assert_success(output: &Output) {
    assert!(output.status.success(), "stdout:\n{}\nstderr:\n{}", String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
}

fn assert_rejected(output: &Output) {
    assert!(!output.status.success(), "unsupported or malformed input was accepted: stdout={} stderr={}", String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
}

fn assert_output_contract(output: &Path, index: &str, expected_record: &str) -> Result<String> {
    assert!(output.exists());
    assert!(output.to_string_lossy().ends_with(".vcf.gz"));
    let index_path = PathBuf::from(format!("{}.{}", output.display(), index));
    assert!(index_path.exists(), "missing {} sidecar", index_path.display());
    let text = read_bgzf(output)?;
    assert!(text.contains("##fileformat=VCFv4.3"));
    assert!(text.contains(expected_record));
    Ok(text)
}

#[test]
fn annotate_matrix_accepts_plain_and_bgzf_inputs_and_writes_both_indexes() -> Result<()> {
    let dir = tempdir()?;
    let database = dir.path().join("database.vcf");
    write_plain(&database, &format!("{HEADER}chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25\n"))?;

    for (suffix, make_input) in [("plain", write_plain as fn(&Path, &str) -> Result<()>), ("bgzf", write_bgzf)] {
        let input = dir.path().join(format!("input-{suffix}.vcf{}", if suffix == "bgzf" { ".bgzf" } else { "" }));
        make_input(&input, &format!("{HEADER}{RECORDS}"))?;
        for index in ["csi", "tbi"] {
            let output = dir.path().join(format!("annotated-{suffix}-{index}.vcf.gz"));
            let output_result = run(&["annotate", "--input", input.to_str().unwrap(), "--database", &format!("db={}", database.display()), "--annotation", "db:AF=db_AF", "--output", output.to_str().unwrap(), "--index-type", index])?;
            assert_success(&output_result);
            let text = assert_output_contract(&output, index, "chr1\t10\trs10\tA\tC")?;
            assert!(text.contains("##INFO=<ID=db_AF"));
            assert!(text.contains("##source=contract-fixture"));
        }
    }
    Ok(())
}

#[test]
fn filter_matrix_preserves_header_and_selected_source_record() -> Result<()> {
    let dir = tempdir()?;
    for (suffix, make_input) in [("plain", write_plain as fn(&Path, &str) -> Result<()>), ("bgzf", write_bgzf)] {
        let input = dir.path().join(format!("input-{suffix}.vcf{}", if suffix == "bgzf" { ".gz" } else { "" }));
        make_input(&input, &format!("{HEADER}{RECORDS}"))?;
        let output = dir.path().join(format!("filtered-{suffix}.vcf.gz"));
        let result = run(&["filter", "--input", input.to_str().unwrap(), "--max-info", "AF=0.10", "--output", output.to_str().unwrap(), "--index-type", "tbi"])?;
        assert_success(&result);
        let text = assert_output_contract(&output, "tbi", "chr1\t10\trs10\tA\tC")?;
        assert!(text.contains("##source=contract-fixture"));
        assert!(!text.contains("chr1\t11\trs11"));
    }
    Ok(())
}

#[test]
fn intersect_matrix_matches_literal_keys_and_preserves_selected_source_header() -> Result<()> {
    let dir = tempdir()?;
    let left = dir.path().join("left.vcf");
    let right = dir.path().join("right.vcf.gz");
    write_plain(&left, &format!("{HEADER}{RECORDS}"))?;
    write_bgzf(&right, &format!("{HEADER}chr1\t10\tother\tA\tC\t.\tPASS\tDB=1\n"))?;

    for index in ["csi", "tbi"] {
        let output = dir.path().join(format!("shared-{index}.vcf.gz"));
        let result = run(&["intersect", "--left", left.to_str().unwrap(), "--right", right.to_str().unwrap(), "--mode", "shared", "--output", output.to_str().unwrap(), "--index-type", index])?;
        assert_success(&result);
        let text = assert_output_contract(&output, index, "chr1\t10\trs10\tA\tC")?;
        assert!(text.contains("##source=contract-fixture"));
        assert!(!text.contains("DB=1"), "intersect must emit the selected source record");
    }
    Ok(())
}

#[test]
fn annotate_reference_normalization_matches_case_insensitive_alleles() -> Result<()> {
    let dir = tempdir()?;
    let reference = write_reference(dir.path())?;
    let input = dir.path().join("lower.vcf");
    let database = dir.path().join("db.vcf");
    let output = dir.path().join("normalized.vcf.gz");
    write_plain(&input, &format!("{HEADER}chr1\t2\t.\ta\tc\t.\tPASS\t.\n"))?;
    write_plain(&database, &format!("{HEADER}chr1\t2\t.\tA\tC\t.\tPASS\tAF=0.75\n"))?;
    let result = run(&["annotate", "--input", input.to_str().unwrap(), "--database", &format!("db={}", database.display()), "--annotation", "db:AF=db_AF", "--reference", reference.to_str().unwrap(), "--output", output.to_str().unwrap()])?;
    assert_success(&result);
    let text = assert_output_contract(&output, "csi", "chr1\t2\t.\ta\tc")?;
    assert!(text.contains("db_AF=0.75"));
    Ok(())
}
#[test]
fn malformed_vcf_and_bcf_are_rejected_by_transforming_commands() -> Result<()> {
    let dir = tempdir()?;
    let malformed = dir.path().join("malformed.vcf");
    let bcf = dir.path().join("input.bcf");
    write_plain(&malformed, &format!("{HEADER}chr1\tbad\n"))?;
    fs::write(&bcf, b"BCF\x02\x02not-a-text-vcf")?;
    let database = dir.path().join("db.vcf");
    write_plain(&database, &format!("{HEADER}chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.1\n"))?;

    let database_arg = format!("db={}", database.display());
    let cases = [
        vec!["annotate", "--database", database_arg.as_str(), "--annotation", "db:AF=x"],
        vec!["filter", "--max-info", "AF=1"],
        vec!["intersect", "--mode", "shared"],
    ];
    for input in [&malformed, &bcf] {
        for case in &cases {
            let output = dir.path().join("rejected.vcf.gz");
            let mut args = case.clone();
            if case[0] == "annotate" { args.extend(["--input", input.to_str().unwrap()]); }
            else if case[0] == "filter" { args.extend(["--input", input.to_str().unwrap()]); }
            else { args.extend(["--left", input.to_str().unwrap(), "--right", input.to_str().unwrap()]); }
            args.extend(["--output", output.to_str().unwrap()]);
            assert_rejected(&run(&args)?);
        }
    }
    Ok(())
}
