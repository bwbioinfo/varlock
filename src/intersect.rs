use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, Write},
    path::Path,
    time::Instant,
};

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::{
    ExecutionContext, IntersectArgs, IntersectMode, log_verbose,
    vcf::{IndexRecord, OutputIndex, VariantKey, open_text_reader, variant_keys_from_fields},
};

#[derive(Default)]
struct IntersectMetrics {
    left_records: usize,
    right_records: usize,
    output_records: usize,
}

pub(crate) fn run(args: IntersectArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    let left_record_count = count_records(&args.left)?;
    let right_record_count = count_records(&args.right)?;
    let right_keys = load_variant_keys(&args.right)?;
    let left_keys = if matches!(args.mode, IntersectMode::RightOnly) {
        load_variant_keys(&args.left)?
    } else {
        HashSet::new()
    };
    let mut metrics = match args.mode {
        IntersectMode::Shared => write_selected_records(
            &args.left,
            &args.output,
            args.index_type,
            &right_keys,
            MatchSense::Present,
            "both",
        )?,
        IntersectMode::LeftOnly => write_selected_records(
            &args.left,
            &args.output,
            args.index_type,
            &right_keys,
            MatchSense::Absent,
            "left",
        )?,
        IntersectMode::RightOnly => write_selected_records(
            &args.right,
            &args.output,
            args.index_type,
            &left_keys,
            MatchSense::Absent,
            "right",
        )?,
    };
    metrics.left_records = left_record_count;
    metrics.right_records = right_record_count;

    log_verbose(
        ctx,
        format!(
            "intersect stage=done mode={:?} left_records={} right_records={} output_records={} output={} elapsed={:.2?}",
            args.mode,
            metrics.left_records,
            metrics.right_records,
            metrics.output_records,
            args.output.display(),
            started.elapsed()
        ),
    );
    Ok(())
}

#[derive(Clone, Copy)]
enum MatchSense {
    Present,
    Absent,
}

fn load_variant_keys(path: &Path) -> Result<HashSet<VariantKey>> {
    let mut reader =
        open_text_reader(path).with_context(|| format!("failed to open VCF {}", path.display()))?;
    let mut keys = HashSet::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            keys.extend(variant_keys_from_fields(&fields));
        }
        line.clear();
    }
    Ok(keys)
}

fn count_records(path: &Path) -> Result<usize> {
    let mut reader =
        open_text_reader(path).with_context(|| format!("failed to open VCF {}", path.display()))?;
    let mut count = 0usize;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            count += 1;
        }
        line.clear();
    }
    Ok(count)
}

fn write_selected_records(
    input: &Path,
    output: &Path,
    index_type: crate::IndexType,
    comparison_keys: &HashSet<VariantKey>,
    sense: MatchSense,
    membership: &str,
) -> Result<IntersectMetrics> {
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let mut reader = open_text_reader(input)
        .with_context(|| format!("failed to open input VCF {}", input.display()))?;
    let output_file = File::create(output)
        .with_context(|| format!("failed to create output {}", output.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    let mut output_index = OutputIndex::new(index_type);
    let mut metrics = IntersectMetrics::default();
    let mut wrote_header = false;
    let mut wrote_varlock_header = false;

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            if !wrote_varlock_header {
                write_varlock_intersect_headers(&mut writer)?;
                wrote_varlock_header = true;
            }
            writeln!(writer, "{trimmed}")?;
            wrote_header = true;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.left_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let keys = variant_keys_from_fields(&fields);
            let is_match = keys.iter().any(|key| comparison_keys.contains(key));
            let keep = match sense {
                MatchSense::Present => is_match,
                MatchSense::Absent => !is_match,
            };
            if keep {
                let annotated = add_membership_info(&fields, membership);
                let index_record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{annotated}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&index_record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }

    if !wrote_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

fn write_varlock_intersect_headers<W: Write>(writer: &mut W) -> Result<()> {
    writeln!(
        writer,
        "##INFO=<ID=VARLOCK_SET,Number=1,Type=String,Description=\"varlock intersect membership: left, right, or both\">"
    )?;
    Ok(())
}

fn add_membership_info(fields: &[&str], membership: &str) -> String {
    let mut out = fields.to_vec();
    let mut info = if out[7] == "." || out[7].is_empty() {
        String::new()
    } else {
        out[7].to_string()
    };
    if !info.is_empty() {
        info.push(';');
    }
    info.push_str("VARLOCK_SET=");
    info.push_str(membership);
    out[7] = &info;
    out.join("\t")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn intersect_shared_writes_left_records_with_membership() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let output = dir.path().join("shared.vcf.gz");
        write_vcf(
            &left,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.2\n",
        )?;
        write_vcf(&right, "chr1\t10\t.\tA\tC\t.\tPASS\t.\n")?;

        let args = IntersectArgs {
            left,
            right,
            mode: IntersectMode::Shared,
            output: output.clone(),
            index_type: crate::IndexType::Csi,
        };
        run(
            args,
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains("##INFO=<ID=VARLOCK_SET"));
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tVARLOCK_SET=both"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(dir.path().join("shared.vcf.gz.csi").exists());
        Ok(())
    }

    #[test]
    fn intersect_left_only_and_right_only_choose_output_side() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let left_only = dir.path().join("left_only.vcf.gz");
        let right_only = dir.path().join("right_only.vcf.gz");
        write_vcf(
            &left,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t11\t.\tA\tG\t.\tPASS\t.\n",
        )?;
        write_vcf(
            &right,
            "chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t12\t.\tA\tT\t.\tPASS\t.\n",
        )?;

        run(
            IntersectArgs {
                left: left.clone(),
                right: right.clone(),
                mode: IntersectMode::LeftOnly,
                output: left_only.clone(),
                index_type: crate::IndexType::Csi,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        run(
            IntersectArgs {
                left,
                right,
                mode: IntersectMode::RightOnly,
                output: right_only.clone(),
                index_type: crate::IndexType::Csi,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let left_text = read_bgzip(&left_only)?;
        assert!(left_text.contains("chr1\t11\t.\tA\tG\t.\tPASS\tVARLOCK_SET=left"));
        assert!(!left_text.contains("chr1\t12\t.\tA\tT"));

        let right_text = read_bgzip(&right_only)?;
        assert!(right_text.contains("chr1\t12\t.\tA\tT\t.\tPASS\tVARLOCK_SET=right"));
        assert!(!right_text.contains("chr1\t11\t.\tA\tG"));
        Ok(())
    }

    #[test]
    fn intersect_matches_any_alt_in_multialt_record() -> Result<()> {
        let dir = tempdir()?;
        let left = dir.path().join("left.vcf");
        let right = dir.path().join("right.vcf");
        let output = dir.path().join("shared.vcf.gz");
        write_vcf(&left, "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\n")?;
        write_vcf(&right, "chr1\t10\t.\tA\tG\t.\tPASS\t.\n")?;

        run(
            IntersectArgs {
                left,
                right,
                mode: IntersectMode::Shared,
                output: output.clone(),
                index_type: crate::IndexType::Csi,
            },
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;

        let text = read_bgzip(&output)?;
        assert!(text.contains("chr1\t10\t.\tA\tC,G\t.\tPASS\tVARLOCK_SET=both"));
        Ok(())
    }

    fn write_vcf(path: &Path, records: &str) -> Result<()> {
        std::fs::write(
            path,
            format!(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"
            ),
        )?;
        Ok(())
    }

    fn read_bgzip(path: &Path) -> Result<String> {
        let mut reader = bgzf::io::Reader::new(File::open(path)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        Ok(text)
    }
}
