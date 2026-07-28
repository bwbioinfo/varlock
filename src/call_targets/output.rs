use std::{
    collections::BTreeMap,
    fs::File,
    io::Write,
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::{
    self,
    index::reference_sequence::bin::Chunk,
    index::reference_sequence::index::BinnedIndex,
    index::{Header as TabixHeader, header::Format as TabixFormat},
};
use noodles_tabix as tabix;

use crate::{CallTargetsArgs, ExecutionContext, IndexType, log_verbose};

use super::reference::{FastaIndex, open_fasta_index};
use super::types::{PairedCallingConfig, SiteCounts, SiteKey, base_index};

pub(crate) struct CallTargetsOutputState {
    writer: bgzf::io::Writer<File>,
    fasta: FastaIndex,
    csi_indexer: Option<binning_index::Indexer<BinnedIndex>>,
    tbi_indexer: Option<tabix::index::Indexer>,
    output_path: PathBuf,
    visited_sites: usize,
    written_variants: usize,
    write_started: Instant,
    write_metrics: WriteMetrics,
}

#[derive(Default)]
struct WriteMetrics {
    fasta_lookup_ns: u64,
    filter_eval_ns: u64,
    vcf_write_ns: u64,
    index_update_ns: u64,
}

pub(crate) fn write_call_targets_output(
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
    prepared_reference: &Path,
    output_path: &Path,
    ref_names: &[String],
    sample_names: &[String],
    paired: Option<&PairedCallingConfig>,
    counts: BTreeMap<SiteKey, SiteCounts>,
) -> Result<()> {
    let total_sites = counts.len();
    let mut state = begin_call_targets_output(
        args,
        ctx,
        label,
        prepared_reference,
        output_path,
        ref_names,
        sample_names,
        paired,
    )?;
    write_call_targets_output_chunk(
        &mut state,
        args,
        ctx,
        label,
        ref_names,
        paired,
        counts,
        Some(total_sites),
    )?;
    finish_call_targets_output(state, args, ctx, label, ref_names)
}

pub(crate) fn begin_call_targets_output(
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
    prepared_reference: &Path,
    output_path: &Path,
    ref_names: &[String],
    sample_names: &[String],
    paired: Option<&PairedCallingConfig>,
) -> Result<CallTargetsOutputState> {
    let stage_started = Instant::now();
    let fasta = open_fasta_index(prepared_reference)?;
    log_verbose(
        ctx,
        format!(
            "{label} stage=open_reference elapsed={:.2?}",
            stage_started.elapsed()
        ),
    );

    if let Some(parent) = output_path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let output_file = File::create(output_path)
        .with_context(|| format!("failed to create output {}", output_path.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    write_vcf_header(&mut writer, args, ref_names, sample_names, &fasta, paired)?;

    let csi_indexer = matches!(args.index_type, IndexType::Csi)
        .then(binning_index::Indexer::<BinnedIndex>::default);
    let tbi_indexer = matches!(args.index_type, IndexType::Tbi).then(|| {
        let mut indexer = tabix::index::Indexer::default();
        let header = TabixHeader::builder()
            .set_format(TabixFormat::Vcf)
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'#')
            .set_line_skip_count(0)
            .build();
        indexer.set_header(header);
        indexer
    });

    Ok(CallTargetsOutputState {
        writer,
        fasta,
        csi_indexer,
        tbi_indexer,
        output_path: output_path.to_path_buf(),
        visited_sites: 0,
        written_variants: 0,
        write_started: Instant::now(),
        write_metrics: WriteMetrics::default(),
    })
}

pub(crate) fn write_call_targets_output_chunk(
    state: &mut CallTargetsOutputState,
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
    ref_names: &[String],
    paired: Option<&PairedCallingConfig>,
    counts: BTreeMap<SiteKey, SiteCounts>,
    total_sites: Option<usize>,
) -> Result<()> {
    for (key, site_counts) in counts {
        state.visited_sites += 1;
        let ref_name = ref_names
            .get(key.reference_sequence_id)
            .context("reference sequence id out of range")?;

        let fasta_lookup_started = Instant::now();
        let ref_base = state.fasta.fetch_base(ref_name, key.position)?;
        state.write_metrics.fasta_lookup_ns = state
            .write_metrics
            .fasta_lookup_ns
            .saturating_add(elapsed_ns(fasta_lookup_started));

        let filter_eval_started = Instant::now();
        if base_index(ref_base).is_none() {
            state.write_metrics.filter_eval_ns = state
                .write_metrics
                .filter_eval_ns
                .saturating_add(elapsed_ns(filter_eval_started));
            continue;
        }
        let (alt_base, total_alt) = choose_alt(ref_base, &site_counts, args.min_alt_count)?;
        if alt_base.is_none() {
            state.write_metrics.filter_eval_ns = state
                .write_metrics
                .filter_eval_ns
                .saturating_add(elapsed_ns(filter_eval_started));
            continue;
        }

        let alt_base = alt_base.expect("checked");
        let total_dp: u32 = site_counts
            .per_sample
            .iter()
            .map(|counts| counts.iter().sum::<u32>())
            .sum();
        if total_dp == 0 {
            state.write_metrics.filter_eval_ns = state
                .write_metrics
                .filter_eval_ns
                .saturating_add(elapsed_ns(filter_eval_started));
            continue;
        }

        let alt_fraction = total_alt as f64 / total_dp as f64;
        if alt_fraction < args.min_alt_fraction {
            state.write_metrics.filter_eval_ns = state
                .write_metrics
                .filter_eval_ns
                .saturating_add(elapsed_ns(filter_eval_started));
            continue;
        }
        let paired_call = match paired {
            Some(config) => match evaluate_paired_call(config, alt_base, &site_counts)? {
                Some(call) => Some(call),
                None => {
                    state.write_metrics.filter_eval_ns = state
                        .write_metrics
                        .filter_eval_ns
                        .saturating_add(elapsed_ns(filter_eval_started));
                    continue;
                }
            },
            None => None,
        };
        state.write_metrics.filter_eval_ns = state
            .write_metrics
            .filter_eval_ns
            .saturating_add(elapsed_ns(filter_eval_started));

        let vcf_write_started = Instant::now();
        let chunk_start = state.writer.virtual_position();
        write_vcf_record(
            &mut state.writer,
            ref_name,
            key.position,
            ref_base,
            alt_base,
            total_dp,
            paired_call.as_ref(),
            &site_counts,
        )?;
        state.written_variants += 1;
        let chunk_end = state.writer.virtual_position();
        state.write_metrics.vcf_write_ns = state
            .write_metrics
            .vcf_write_ns
            .saturating_add(elapsed_ns(vcf_write_started));

        let index_update_started = Instant::now();
        let chunk = Chunk::new(chunk_start, chunk_end);
        let start = Position::try_from(key.position as usize)
            .context("invalid VCF position for indexing")?;
        if let Some(indexer) = state.csi_indexer.as_mut() {
            indexer
                .add_record(Some((key.reference_sequence_id, start, start, true)), chunk)
                .context("failed to update CSI index")?;
        }
        if let Some(indexer) = state.tbi_indexer.as_mut() {
            indexer
                .add_record(ref_name, start, start, chunk)
                .context("failed to update TBI index")?;
        }
        state.write_metrics.index_update_ns = state
            .write_metrics
            .index_update_ns
            .saturating_add(elapsed_ns(index_update_started));

        let reached_total = total_sites
            .map(|expected| state.visited_sites == expected)
            .unwrap_or(false);
        if ctx.verbose > 0 && (state.visited_sites.is_multiple_of(100_000) || reached_total) {
            eprintln!(
                "[{label}] write progress: sites={}/{} variants={} elapsed={:.2?}",
                state.visited_sites,
                total_sites.unwrap_or(state.visited_sites),
                state.written_variants,
                state.write_started.elapsed()
            );
        }
    }

    Ok(())
}

pub(crate) fn finish_call_targets_output(
    mut state: CallTargetsOutputState,
    args: &CallTargetsArgs,
    ctx: &ExecutionContext,
    label: &str,
    ref_names: &[String],
) -> Result<()> {
    let writer_finish_started = Instant::now();
    state
        .writer
        .try_finish()
        .context("failed to finish writing VCF")?;
    if ctx.verbose > 0 {
        eprintln!(
            "[{label}] writer finalize: bgzf_try_finish_elapsed={:.2?}",
            writer_finish_started.elapsed()
        );
    }
    log_verbose(
        ctx,
        format!(
            "{label} stage=write_vcf variants={} elapsed={:.2?}",
            state.written_variants,
            state.write_started.elapsed()
        ),
    );
    if ctx.verbose > 0 {
        eprintln!(
            "[{label}] write breakdown: fasta_lookup_ms={:.2} filter_eval_ms={:.2} vcf_write_ms={:.2} index_update_ms={:.2}",
            state.write_metrics.fasta_lookup_ns as f64 / 1_000_000.0,
            state.write_metrics.filter_eval_ns as f64 / 1_000_000.0,
            state.write_metrics.vcf_write_ns as f64 / 1_000_000.0,
            state.write_metrics.index_update_ns as f64 / 1_000_000.0,
        );
    }

    match args.index_type {
        IndexType::Csi => {
            let index_build_started = Instant::now();
            let indexer = state.csi_indexer.context("CSI indexer missing")?;
            let index = indexer.build(ref_names.len());
            if ctx.verbose > 0 {
                eprintln!(
                    "[{label}] index finalize: type=csi build_elapsed={:.2?}",
                    index_build_started.elapsed()
                );
            }
            let index_path = state
                .output_path
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name| state.output_path.with_file_name(format!("{}.csi", name)))
                .context("invalid output path")?;
            let index_write_started = Instant::now();
            let index_file = File::create(&index_path)
                .with_context(|| format!("failed to create index {}", index_path.display()))?;
            let mut csi_writer = noodles_csi::io::Writer::new(index_file);
            csi_writer
                .write_index(&index)
                .context("failed to write CSI index")?;
            if ctx.verbose > 0 {
                eprintln!(
                    "[{label}] index finalize: type=csi write_elapsed={:.2?}",
                    index_write_started.elapsed()
                );
            }
        }
        IndexType::Tbi => {
            let index_build_started = Instant::now();
            let indexer = state.tbi_indexer.context("TBI indexer missing")?;
            let index = indexer.build();
            if ctx.verbose > 0 {
                eprintln!(
                    "[{label}] index finalize: type=tbi build_elapsed={:.2?}",
                    index_build_started.elapsed()
                );
            }
            let index_path = state
                .output_path
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name| state.output_path.with_file_name(format!("{}.tbi", name)))
                .context("invalid output path")?;
            let index_write_started = Instant::now();
            let index_file = File::create(&index_path)
                .with_context(|| format!("failed to create index {}", index_path.display()))?;
            let mut tbi_writer = tabix::io::Writer::new(index_file);
            tbi_writer
                .write_index(&index)
                .context("failed to write TBI index")?;
            if ctx.verbose > 0 {
                eprintln!(
                    "[{label}] index finalize: type=tbi write_elapsed={:.2?}",
                    index_write_started.elapsed()
                );
            }
        }
    }

    Ok(())
}

fn write_vcf_header<W: Write>(
    writer: &mut bgzf::io::Writer<W>,
    args: &CallTargetsArgs,
    ref_names: &[String],
    sample_names: &[String],
    fasta: &FastaIndex,
    paired: Option<&PairedCallingConfig>,
) -> Result<()> {
    let reference = args
        .reference
        .as_ref()
        .context("--reference is required for call-targets output")?;
    writeln!(writer, "##fileformat=VCFv4.3")?;
    writeln!(writer, "##source=varlock call-targets")?;
    writeln!(writer, "##reference={}", reference.display())?;
    for (name, length) in ref_names.iter().zip(fasta.reference_lengths(ref_names)?) {
        writeln!(writer, "##contig=<ID={},length={}>", name, length)?;
    }
    writeln!(
        writer,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )?;
    if paired.is_some() {
        writeln!(
            writer,
            "##INFO=<ID=PAIR,Number=1,Type=String,Description=\"Paired tumor/normal sample names\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Variant passed paired tumor/normal filters\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=TUMOR_AF,Number=1,Type=Float,Description=\"Tumor alternate allele fraction\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=NORMAL_AF,Number=1,Type=Float,Description=\"Normal alternate allele fraction\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=TUMOR_ALT_COUNT,Number=1,Type=Integer,Description=\"Tumor alternate allele count\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=NORMAL_ALT_COUNT,Number=1,Type=Integer,Description=\"Normal alternate allele count\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=TUMOR_DP,Number=1,Type=Integer,Description=\"Tumor read depth\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=NORMAL_DP,Number=1,Type=Integer,Description=\"Normal read depth\">"
        )?;
    }
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">"
    )?;
    write!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    )?;
    for sample in sample_names {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;
    Ok(())
}

fn write_vcf_record<W: Write>(
    writer: &mut bgzf::io::Writer<W>,
    chrom: &str,
    pos: u32,
    ref_base: u8,
    alt_base: u8,
    total_dp: u32,
    paired_call: Option<&PairedCallInfo>,
    site_counts: &SiteCounts,
) -> Result<()> {
    write!(
        writer,
        "{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={}",
        chrom, pos, ref_base as char, alt_base as char, total_dp
    )?;
    if let Some(call) = paired_call {
        write!(
            writer,
            ";PAIR={};SOMATIC;TUMOR_AF={:.6};NORMAL_AF={:.6};TUMOR_ALT_COUNT={};NORMAL_ALT_COUNT={};TUMOR_DP={};NORMAL_DP={}",
            call.pair,
            call.tumor_af,
            call.normal_af,
            call.tumor_alt_count,
            call.normal_alt_count,
            call.tumor_dp,
            call.normal_dp
        )?;
    }
    write!(writer, "\tGT:DP:AD")?;
    let ref_idx = base_index(ref_base).context("invalid reference base")?;
    let alt_idx = base_index(alt_base).context("invalid alt base")?;

    for sample in &site_counts.per_sample {
        let dp = sample.iter().sum::<u32>();
        let ref_count = sample[ref_idx];
        let alt_count = sample[alt_idx];
        let gt = if dp == 0 {
            "./.".to_string()
        } else if alt_count == 0 {
            "0/0".to_string()
        } else if ref_count == 0 {
            "1/1".to_string()
        } else {
            "0/1".to_string()
        };
        write!(writer, "\t{}:{}:{},{}", gt, dp, ref_count, alt_count)?;
    }
    writeln!(writer)?;
    Ok(())
}

#[derive(Clone, Debug)]
struct PairedCallInfo {
    pair: String,
    tumor_dp: u32,
    tumor_alt_count: u32,
    tumor_af: f64,
    normal_dp: u32,
    normal_alt_count: u32,
    normal_af: f64,
}

fn evaluate_paired_call(
    config: &PairedCallingConfig,
    alt_base: u8,
    site_counts: &SiteCounts,
) -> Result<Option<PairedCallInfo>> {
    let alt_idx = base_index(alt_base).context("invalid paired alt base")?;
    let tumor = site_counts
        .per_sample
        .get(config.tumor_index)
        .context("paired tumor sample index out of range")?;
    let normal = site_counts
        .per_sample
        .get(config.normal_index)
        .context("paired normal sample index out of range")?;

    let tumor_dp = tumor.iter().sum::<u32>();
    let normal_dp = normal.iter().sum::<u32>();
    if config
        .normal_min_depth
        .is_some_and(|min_depth| normal_dp < min_depth)
    {
        return Ok(None);
    }

    let tumor_alt_count = tumor[alt_idx];
    let normal_alt_count = normal[alt_idx];
    let tumor_af = allele_fraction(tumor_alt_count, tumor_dp);
    let normal_af = allele_fraction(normal_alt_count, normal_dp);

    if tumor_alt_count < config.tumor_min_alt_count {
        return Ok(None);
    }
    if tumor_af < config.tumor_min_alt_fraction {
        return Ok(None);
    }
    if normal_alt_count > config.normal_max_alt_count {
        return Ok(None);
    }
    if normal_af > config.normal_max_alt_fraction {
        return Ok(None);
    }

    Ok(Some(PairedCallInfo {
        pair: format!("{}|{}", config.roles.tumor, config.roles.normal),
        tumor_dp,
        tumor_alt_count,
        tumor_af,
        normal_dp,
        normal_alt_count,
        normal_af,
    }))
}

fn allele_fraction(alt_count: u32, depth: u32) -> f64 {
    if depth == 0 {
        0.0
    } else {
        alt_count as f64 / depth as f64
    }
}

fn choose_alt(ref_base: u8, counts: &SiteCounts, min_alt_count: u32) -> Result<(Option<u8>, u32)> {
    let mut totals = [0u32; 4];
    for sample in &counts.per_sample {
        for (i, value) in sample.iter().enumerate() {
            totals[i] += value;
        }
    }

    let Some(ref_idx) = base_index(ref_base) else {
        return Ok((None, 0));
    };
    let mut best_idx = None;
    let mut best_count = 0u32;
    for (i, count) in totals.iter().enumerate() {
        if i == ref_idx {
            continue;
        }
        if *count > best_count {
            best_count = *count;
            best_idx = Some(i);
        }
    }

    if best_count < min_alt_count {
        return Ok((None, 0));
    }

    let alt_base = match best_idx {
        Some(0) => b'A',
        Some(1) => b'C',
        Some(2) => b'G',
        Some(3) => b'T',
        _ => return Ok((None, 0)),
    };

    Ok((Some(alt_base), best_count))
}

fn elapsed_ns(start: Instant) -> u64 {
    start.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "wgpu")]
    use super::super::types::SiteKey;
    use super::super::types::{PairedCallingConfig, PairedSampleRoles, SiteCounts};
    #[cfg(feature = "wgpu")]
    use super::write_call_targets_output;
    use super::{choose_alt, evaluate_paired_call};
    #[cfg(feature = "wgpu")]
    use crate::{CallTargetsArgs, ExecutionContext, IndexType};
    #[cfg(feature = "wgpu")]
    use anyhow::Result;
    #[cfg(feature = "wgpu")]
    use noodles_bgzf as bgzf;
    #[cfg(feature = "wgpu")]
    use std::{collections::BTreeMap, fs::File, io::Read};
    #[cfg(feature = "wgpu")]
    use tempfile::tempdir;

    fn counts(per_sample: Vec<[u32; 4]>) -> SiteCounts {
        SiteCounts { per_sample }
    }

    fn paired_config() -> PairedCallingConfig {
        PairedCallingConfig {
            roles: PairedSampleRoles {
                tumor: "tumor".to_string(),
                normal: "normal".to_string(),
            },
            tumor_index: 0,
            normal_index: 1,
            tumor_min_alt_count: 3,
            tumor_min_alt_fraction: 0.2,
            normal_max_alt_count: 1,
            normal_max_alt_fraction: 0.05,
            normal_min_depth: Some(10),
        }
    }

    #[cfg(feature = "wgpu")]
    fn call_targets_args(reference: std::path::PathBuf) -> CallTargetsArgs {
        CallTargetsArgs {
            inputs: Vec::new(),
            bamlist: None,
            reference: Some(reference),
            targets: None,
            output: None,
            rg_map: None,
            split_by: crate::SplitBy::Sm,
            #[cfg(feature = "wgpu")]
            gpu: crate::GpuArgs::default(),
            index_type: IndexType::Csi,
            min_mapq: 20,
            min_baseq: 20,
            min_alt_count: 1,
            min_alt_fraction: 0.0,
            pair: Some("tumor=tumor,normal=normal".to_string()),
            tumor_min_alt_count: 3,
            tumor_min_alt_fraction: 0.2,
            normal_max_alt_count: 1,
            normal_max_alt_fraction: 0.05,
            normal_min_depth: Some(10),
            max_depth: 100_000,
        }
    }

    #[test]
    fn choose_alt_selects_highest_count_non_ref() {
        // ref=A(0): C=3, G=1, T=0 → best alt is C
        let c = counts(vec![[10, 3, 1, 0]]);
        let (alt, count) = choose_alt(b'A', &c, 1).unwrap();
        assert_eq!(alt, Some(b'C'));
        assert_eq!(count, 3);
    }

    #[test]
    fn choose_alt_returns_none_below_min_count() {
        let c = counts(vec![[10, 2, 0, 0]]);
        let (alt, _) = choose_alt(b'A', &c, 5).unwrap();
        assert_eq!(alt, None);
    }

    #[test]
    fn choose_alt_sums_counts_across_samples() {
        // Two samples each contributing 3 to C → total 6
        let c = counts(vec![[5, 3, 0, 0], [5, 3, 0, 0]]);
        let (alt, count) = choose_alt(b'A', &c, 1).unwrap();
        assert_eq!(alt, Some(b'C'));
        assert_eq!(count, 6);
    }

    #[test]
    fn choose_alt_returns_none_for_non_acgt_ref() {
        let c = counts(vec![[0, 5, 0, 0]]);
        let (alt, count) = choose_alt(b'N', &c, 1).unwrap();
        assert_eq!(alt, None);
        assert_eq!(count, 0);
    }

    #[test]
    fn choose_alt_returns_none_when_only_ref_present() {
        let c = counts(vec![[10, 0, 0, 0]]);
        let (alt, _) = choose_alt(b'A', &c, 1).unwrap();
        assert_eq!(alt, None);
    }

    #[test]
    fn choose_alt_respects_ref_base_when_selecting_alt() {
        // ref=C(1), A=5 is highest non-ref → alt should be A
        let c = counts(vec![[5, 10, 1, 0]]);
        let (alt, count) = choose_alt(b'C', &c, 1).unwrap();
        assert_eq!(alt, Some(b'A'));
        assert_eq!(count, 5);
    }

    #[test]
    fn evaluate_paired_call_accepts_tumor_supported_normal_clean_alt() {
        let c = counts(vec![[10, 4, 0, 0], [20, 1, 0, 0]]);
        let call = evaluate_paired_call(&paired_config(), b'C', &c)
            .unwrap()
            .unwrap();

        assert_eq!(call.tumor_dp, 14);
        assert_eq!(call.tumor_alt_count, 4);
        assert!((call.tumor_af - 0.285714).abs() < 0.000001);
        assert_eq!(call.normal_dp, 21);
        assert_eq!(call.normal_alt_count, 1);
        assert!((call.normal_af - 0.047619).abs() < 0.000001);
    }

    #[test]
    fn evaluate_paired_call_rejects_weak_tumor_or_contaminated_normal() {
        let weak_tumor = counts(vec![[10, 2, 0, 0], [20, 0, 0, 0]]);
        assert!(
            evaluate_paired_call(&paired_config(), b'C', &weak_tumor)
                .unwrap()
                .is_none()
        );

        let contaminated_normal = counts(vec![[10, 4, 0, 0], [20, 2, 0, 0]]);
        assert!(
            evaluate_paired_call(&paired_config(), b'C', &contaminated_normal)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn evaluate_paired_call_rejects_low_normal_depth_when_requested() {
        let c = counts(vec![[10, 4, 0, 0], [5, 0, 0, 0]]);
        assert!(
            evaluate_paired_call(&paired_config(), b'C', &c)
                .unwrap()
                .is_none()
        );
    }

    #[cfg(feature = "wgpu")]
    #[test]
    fn gpu_output_path_writes_paired_headers_and_info() -> Result<()> {
        let dir = tempdir()?;
        let reference = dir.path().join("ref.fa");
        let output = dir.path().join("gpu.paired.vcf.gz");
        std::fs::write(&reference, ">chr1\nACGTACGT\n")?;
        std::fs::write(reference.with_extension("fa.fai"), "chr1\t8\t6\t8\t9\n")?;

        let mut site_counts = BTreeMap::new();
        site_counts.insert(
            SiteKey {
                reference_sequence_id: 0,
                position: 2,
            },
            counts(vec![[0, 10, 4, 0], [0, 20, 1, 0]]),
        );

        let args = call_targets_args(reference.clone());
        let ctx = ExecutionContext {
            verbose: 0,
            threads: 1,
        };
        write_call_targets_output(
            &args,
            &ctx,
            "call_targets_gpu",
            &reference,
            &output,
            &["chr1".to_string()],
            &["tumor".to_string(), "normal".to_string()],
            Some(&paired_config()),
            site_counts,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;

        assert!(text.contains("##INFO=<ID=SOMATIC"));
        assert!(text.contains("##INFO=<ID=TUMOR_AF"));
        assert!(text.contains(
            "chr1\t2\t.\tC\tG\t.\tPASS\tDP=35;PAIR=tumor|normal;SOMATIC;TUMOR_AF=0.285714;NORMAL_AF=0.047619;TUMOR_ALT_COUNT=4;NORMAL_ALT_COUNT=1;TUMOR_DP=14;NORMAL_DP=21\tGT:DP:AD\t0/1:14:10,4\t0/1:21:20,1"
        ));
        assert!(dir.path().join("gpu.paired.vcf.gz.csi").exists());
        Ok(())
    }
}
