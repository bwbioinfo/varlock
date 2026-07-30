use std::{
    cmp::Ordering,
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
use super::types::{
    IndelAllele, IndelCounts, IndelKey, PairedCallingConfig, SiteCounts, SiteKey, base_index,
};

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

pub(crate) struct CallTargetsOutputContext<'a> {
    pub(crate) args: &'a CallTargetsArgs,
    pub(crate) ctx: &'a ExecutionContext,
    pub(crate) label: &'a str,
    pub(crate) ref_names: &'a [String],
    pub(crate) paired: Option<&'a PairedCallingConfig>,
}

pub(crate) fn write_call_targets_output(
    output_context: CallTargetsOutputContext<'_>,
    prepared_reference: &Path,
    output_path: &Path,
    sample_names: &[String],
    counts: BTreeMap<SiteKey, SiteCounts>,
    indel_counts: BTreeMap<IndelKey, IndelCounts>,
) -> Result<()> {
    let mut calls = Vec::with_capacity(counts.len() + indel_counts.len());
    if output_context.args.emit_snvs() {
        calls.extend(
            counts
                .into_iter()
                .map(|(key, counts)| OutputCall::Snv { key, counts }),
        );
    }
    if output_context.args.emit_indels() {
        calls.extend(
            indel_counts
                .into_iter()
                .map(|(key, counts)| OutputCall::Indel { key, counts }),
        );
    }
    calls.sort_unstable_by(OutputCall::compare);
    let total_calls = calls.len();
    let mut state = begin_call_targets_output(
        &output_context,
        prepared_reference,
        output_path,
        sample_names,
    )?;
    for call in calls {
        write_call_targets_output_call(&mut state, &output_context, call, Some(total_calls))?;
    }
    finish_call_targets_output(
        state,
        output_context.args,
        output_context.ctx,
        output_context.label,
        output_context.ref_names,
    )
}

fn begin_call_targets_output(
    output_context: &CallTargetsOutputContext<'_>,
    prepared_reference: &Path,
    output_path: &Path,
    sample_names: &[String],
) -> Result<CallTargetsOutputState> {
    let stage_started = Instant::now();
    let fasta = open_fasta_index(prepared_reference)?;
    log_verbose(
        output_context.ctx,
        format!(
            "{} stage=open_reference elapsed={:.2?}",
            output_context.label,
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
    write_vcf_header(
        &mut writer,
        output_context.args,
        output_context.ref_names,
        sample_names,
        &fasta,
        output_context.paired,
    )?;

    let csi_indexer = matches!(output_context.args.index_type, IndexType::Csi)
        .then(binning_index::Indexer::<BinnedIndex>::default);
    let tbi_indexer = matches!(output_context.args.index_type, IndexType::Tbi).then(|| {
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

enum OutputCall {
    Snv { key: SiteKey, counts: SiteCounts },
    Indel { key: IndelKey, counts: IndelCounts },
}

impl OutputCall {
    fn sort_key(&self) -> (usize, u32, u8) {
        match self {
            Self::Snv { key, .. } => (key.reference_sequence_id, key.position, 0),
            Self::Indel { key, .. } => (key.reference_sequence_id, key.position, 1),
        }
    }

    fn compare(&self, other: &Self) -> Ordering {
        let order = self.sort_key().cmp(&other.sort_key());
        if order != Ordering::Equal {
            return order;
        }
        match (self, other) {
            (Self::Indel { key: left, .. }, Self::Indel { key: right, .. }) => left.cmp(right),
            _ => Ordering::Equal,
        }
    }
}

#[derive(Clone, Copy)]
struct SampleAlleleDepth {
    depth: u32,
    ref_count: u32,
    alt_count: u32,
}

fn write_call_targets_output_call(
    state: &mut CallTargetsOutputState,
    output_context: &CallTargetsOutputContext<'_>,
    call: OutputCall,
    total_sites: Option<usize>,
) -> Result<()> {
    state.visited_sites += 1;
    let (reference_sequence_id, position, ref_bases, alt_bases, sample_depths) = match call {
        OutputCall::Snv { key, counts } => {
            let ref_name = output_context
                .ref_names
                .get(key.reference_sequence_id)
                .context("reference sequence id out of range")?;
            let fasta_lookup_started = Instant::now();
            let ref_base = state.fasta.fetch_base(ref_name, key.position)?;
            state.write_metrics.fasta_lookup_ns = state
                .write_metrics
                .fasta_lookup_ns
                .saturating_add(elapsed_ns(fasta_lookup_started));
            let (Some(alt_base), _) =
                choose_alt(ref_base, &counts, output_context.args.min_alt_count)?
            else {
                return Ok(());
            };
            let ref_idx = base_index(ref_base).context("invalid reference base")?;
            let alt_idx = base_index(alt_base).context("invalid alternate base")?;
            let sample_depths: Vec<SampleAlleleDepth> = counts
                .per_sample
                .iter()
                .map(|sample| SampleAlleleDepth {
                    depth: sample.iter().sum(),
                    ref_count: sample[ref_idx],
                    alt_count: sample[alt_idx],
                })
                .collect();
            (
                key.reference_sequence_id,
                key.position,
                vec![ref_base],
                vec![alt_base],
                sample_depths,
            )
        }
        OutputCall::Indel { key, counts } => {
            let ref_name = output_context
                .ref_names
                .get(key.reference_sequence_id)
                .context("reference sequence id out of range")?;
            let reference_len = match &key.allele {
                IndelAllele::Insertion(_) => 1,
                IndelAllele::Deletion(length) => length
                    .checked_add(1)
                    .context("indel deletion length exceeds VCF coordinate range")?,
            };
            let fasta_lookup_started = Instant::now();
            let ref_bases = state
                .fasta
                .fetch_bases(ref_name, key.position, reference_len)?;
            state.write_metrics.fasta_lookup_ns = state
                .write_metrics
                .fasta_lookup_ns
                .saturating_add(elapsed_ns(fasta_lookup_started));
            if !ref_bases.iter().all(|&base| base_index(base).is_some()) {
                return Ok(());
            }
            let alt_bases = match &key.allele {
                IndelAllele::Insertion(inserted) => {
                    let mut alt = vec![ref_bases[0]];
                    alt.extend(inserted);
                    alt
                }
                IndelAllele::Deletion(_) => vec![ref_bases[0]],
            };
            let sample_depths: Vec<SampleAlleleDepth> = counts
                .per_sample
                .iter()
                .map(|sample| SampleAlleleDepth {
                    depth: sample.iter().sum(),
                    ref_count: sample[0],
                    alt_count: sample[1],
                })
                .collect();
            (
                key.reference_sequence_id,
                key.position,
                ref_bases,
                alt_bases,
                sample_depths,
            )
        }
    };
    let ref_name = output_context
        .ref_names
        .get(reference_sequence_id)
        .context("reference sequence id out of range")?;
    let total_alt = sample_depths
        .iter()
        .map(|sample| sample.alt_count)
        .sum::<u32>();
    let total_dp = sample_depths.iter().map(|sample| sample.depth).sum::<u32>();
    let filter_eval_started = Instant::now();
    if total_alt < output_context.args.min_alt_count
        || total_dp == 0
        || allele_fraction(total_alt, total_dp) < output_context.args.min_alt_fraction
    {
        state.write_metrics.filter_eval_ns = state
            .write_metrics
            .filter_eval_ns
            .saturating_add(elapsed_ns(filter_eval_started));
        return Ok(());
    }
    let paired_call = match output_context.paired {
        Some(config) => match evaluate_paired_depths(config, &sample_depths)? {
            Some(call) => Some(call),
            None => {
                state.write_metrics.filter_eval_ns = state
                    .write_metrics
                    .filter_eval_ns
                    .saturating_add(elapsed_ns(filter_eval_started));
                return Ok(());
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
        VcfRecord {
            chrom: ref_name,
            pos: position,
            ref_bases: &ref_bases,
            alt_bases: &alt_bases,
            total_dp,
            paired_call: paired_call.as_ref(),
            sample_depths: &sample_depths,
        },
    )?;
    state.written_variants += 1;
    let chunk_end = state.writer.virtual_position();
    state.write_metrics.vcf_write_ns = state
        .write_metrics
        .vcf_write_ns
        .saturating_add(elapsed_ns(vcf_write_started));

    let index_update_started = Instant::now();
    let chunk = Chunk::new(chunk_start, chunk_end);
    let start =
        Position::try_from(position as usize).context("invalid VCF position for indexing")?;
    let end_pos = position
        .checked_add(u32::try_from(ref_bases.len()).context("VCF REF length exceeds u32")? - 1)
        .context("VCF reference span exceeds u32")?;
    let end =
        Position::try_from(end_pos as usize).context("invalid VCF end position for indexing")?;
    if let Some(indexer) = state.csi_indexer.as_mut() {
        indexer
            .add_record(Some((reference_sequence_id, start, end, true)), chunk)
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
    if output_context.ctx.verbose > 0
        && (state.visited_sites.is_multiple_of(100_000) || reached_total)
    {
        eprintln!(
            "[{}] write progress: sites={}/{} variants={} elapsed={:.2?}",
            output_context.label,
            state.visited_sites,
            total_sites.unwrap_or(state.visited_sites),
            state.written_variants,
            state.write_started.elapsed()
        );
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

struct VcfRecord<'a> {
    chrom: &'a str,
    pos: u32,
    ref_bases: &'a [u8],
    alt_bases: &'a [u8],
    total_dp: u32,
    paired_call: Option<&'a PairedCallInfo>,
    sample_depths: &'a [SampleAlleleDepth],
}

fn write_vcf_record<W: Write>(
    writer: &mut bgzf::io::Writer<W>,
    record: VcfRecord<'_>,
) -> Result<()> {
    let ref_bases =
        std::str::from_utf8(record.ref_bases).context("invalid VCF reference allele")?;
    let alt_bases =
        std::str::from_utf8(record.alt_bases).context("invalid VCF alternate allele")?;
    write!(
        writer,
        "{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={}",
        record.chrom, record.pos, ref_bases, alt_bases, record.total_dp
    )?;
    if let Some(call) = record.paired_call {
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

    for sample in record.sample_depths {
        let gt = if sample.depth == 0 {
            "./.".to_string()
        } else if sample.alt_count == 0 {
            "0/0".to_string()
        } else if sample.ref_count == 0 {
            "1/1".to_string()
        } else {
            "0/1".to_string()
        };
        write!(
            writer,
            "\t{}:{}:{},{}",
            gt, sample.depth, sample.ref_count, sample.alt_count
        )?;
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

#[cfg_attr(not(test), allow(dead_code))]
fn evaluate_paired_call(
    config: &PairedCallingConfig,
    alt_base: u8,
    site_counts: &SiteCounts,
) -> Result<Option<PairedCallInfo>> {
    let alt_idx = base_index(alt_base).context("invalid paired alt base")?;
    let sample_depths = site_counts
        .per_sample
        .iter()
        .map(|sample| SampleAlleleDepth {
            depth: sample.iter().sum(),
            ref_count: 0,
            alt_count: sample[alt_idx],
        })
        .collect::<Vec<_>>();
    evaluate_paired_depths(config, &sample_depths)
}

fn evaluate_paired_depths(
    config: &PairedCallingConfig,
    sample_depths: &[SampleAlleleDepth],
) -> Result<Option<PairedCallInfo>> {
    let tumor = sample_depths
        .get(config.tumor_index)
        .context("paired tumor sample index out of range")?;
    let normal = sample_depths
        .get(config.normal_index)
        .context("paired normal sample index out of range")?;

    let tumor_dp = tumor.depth;
    let normal_dp = normal.depth;
    if config
        .normal_min_depth
        .is_some_and(|min_depth| normal_dp < min_depth)
    {
        return Ok(None);
    }

    let tumor_alt_count = tumor.alt_count;
    let normal_alt_count = normal.alt_count;
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
    use super::CallTargetsOutputContext;
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
            no_indels: false,
            indels_only: false,
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
            CallTargetsOutputContext {
                args: &args,
                ctx: &ctx,
                label: "call_targets_gpu",
                ref_names: &["chr1".to_string()],
                paired: Some(&paired_config()),
            },
            &reference,
            &output,
            &["tumor".to_string(), "normal".to_string()],
            site_counts,
            BTreeMap::new(),
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
