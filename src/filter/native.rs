//! Schema-aware filter path when either endpoint is BCF.
//!
//! The VCF rendering is only a view for existing predicates. Selected records
//! are passed directly to noodles writers, never reparsed from this view.
use std::{collections::HashMap, fs::File, io::Write as _, path::Path};

use anyhow::{Context, Result, bail};
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_csi::{
    self as csi,
    binning_index::{
        Indexer,
        index::reference_sequence::{bin::Chunk, index::BinnedIndex},
    },
};
use noodles_vcf::{
    self as vcf,
    variant::{Record, io::Write as _},
};

use super::{FilterMetrics, FilterSpec};
use crate::{
    FilterArgs, IndexType, VariantFormat,
    vcf::{
        IndexRecord, OutputIndex, VcfRecord, csi_index_path, open_text_reader, tabix_index_path,
    },
};

pub(super) fn filter(args: &FilterArgs, spec: &FilterSpec) -> Result<FilterMetrics> {
    // Reject incompatible options before opening or creating any output.
    if args.output_format == VariantFormat::Bcf && matches!(args.index_type, IndexType::Tbi) {
        bail!("BCF output requires --index-type csi; TBI is a text VCF index");
    }
    // Never truncate an input alias or leave an old sidecar attached to new data.
    for path in [
        &args.output,
        &csi_index_path(&args.output),
        &tabix_index_path(&args.output),
    ] {
        // symlink_metadata also sees dangling symlinks: these are not new paths.
        match std::fs::symlink_metadata(path) {
            Ok(_) => bail!(
                "native filter output and sidecars must be new paths: {}",
                path.display()
            ),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
    }

    let mut reader: Box<dyn vcf::variant::io::Read> = match args.input_format {
        VariantFormat::Vcf => Box::new(vcf::io::Reader::new(open_text_reader(&args.input)?)),
        // noodles-bcf 0.83 `new` wraps the File in a BGZF reader.
        VariantFormat::Bcf => Box::new(bcf::io::Reader::new(File::open(&args.input)?)),
    };
    let header = reader.read_variant_header().with_context(|| {
        format!(
            "failed to read {:?} header from {}",
            args.input_format,
            args.input.display()
        )
    })?;
    let samples = header.sample_names().iter().cloned().collect::<Vec<_>>();
    spec.validate_samples(&samples)?;
    let sample_index: HashMap<_, _> = samples
        .into_iter()
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();
    let mut output = Output::create(args, &header)?;
    let mut view = vcf::io::Writer::new(Vec::new());
    let mut metrics = FilterMetrics::default();
    for result in reader.variant_records(&header) {
        metrics.input_records += 1;
        let record =
            result.with_context(|| format!("invalid input record {}", metrics.input_records))?;
        view.get_mut().clear();
        view.write_variant_record(&header, record.as_ref())
            .with_context(|| {
                format!("failed to decode record {} for predicates", metrics.input_records)
            })?;
        let line = std::str::from_utf8(view.get_ref())?.trim_end_matches(['\r', '\n']);
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 8 {
            bail!("decoded variant has fewer than eight fields");
        }
        if spec.keep_record(&VcfRecord::new(&fields, &sample_index))? {
            output
                .write(&header, record.as_ref(), line)
                .with_context(|| {
                    format!("failed to write/index selected record {}", metrics.input_records)
                })?;
            metrics.output_records += 1;
        }
    }
    output
        .finish(&args.output)
        .context("failed to finish native filter output/index")?;
    Ok(metrics)
}

enum Output {
    Vcf(bgzf::io::Writer<File>, OutputIndex),
    Bcf {
        writer: bcf::io::Writer<bgzf::io::Writer<File>>,
        indexer: Indexer<BinnedIndex>,
        maps: vcf::header::StringMaps,
        reference_count: usize,
        previous: Option<(usize, noodles_core::Position)>,
    },
}

impl Output {
    fn create(args: &FilterArgs, header: &vcf::Header) -> Result<Self> {
        // Build and validate output dictionaries before creating the data file.
        let maps = vcf::header::StringMaps::try_from(header)?;
        let reference_count = header
            .contigs()
            .keys()
            .filter_map(|name| maps.contigs().get_index_of(name))
            .max()
            .map_or(0, |i| i + 1);
        if let Some(parent) = args.output.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent)?;
        }
        let file = File::options()
            .write(true)
            .create_new(true)
            .open(&args.output)?;
        match args.output_format {
            VariantFormat::Vcf => {
                let mut writer = vcf::io::Writer::new(bgzf::io::Writer::new(file));
                writer.write_header(header)?;
                Ok(Self::Vcf(
                    writer.into_inner(),
                    OutputIndex::new(args.index_type),
                ))
            }
            VariantFormat::Bcf => {
                // `new` provides BGZF transport; `from` would use a raw writer.
                let mut writer = bcf::io::Writer::new(file);
                writer.write_header(header)?;
                Ok(Self::Bcf {
                    writer,
                    // Native BCF CSI has no Tabix text header. Depth 6 covers
                    // the BCF coordinate range, including contigs above 512 Mb.
                    indexer: Indexer::new(14, 6),
                    maps,
                    reference_count,
                    previous: None,
                })
            }
        }
    }

    fn write(&mut self, header: &vcf::Header, record: &dyn Record, line: &str) -> Result<()> {
        match self {
            Self::Vcf(writer, index) => {
                let index_record = IndexRecord::from_line(line)?;
                let start = writer.virtual_position();
                writeln!(writer, "{line}")?;
                index.add_record(&index_record, Chunk::new(start, writer.virtual_position()))?;
            }
            Self::Bcf {
                writer, indexer, maps, previous, ..
            } => {
                let name = record.reference_sequence_name(header)?;
                let id = maps
                    .contigs()
                    .get_index_of(name)
                    .with_context(|| format!("BCF output requires a contig declaration for {name}"))?;
                let start = record
                    .variant_start()
                    .transpose()?
                    .context("BCF indexing requires a POS")?;
                let end = record.variant_end(header)?;
                if end < start {
                    bail!("BCF variant end precedes its start");
                }
                if previous.is_some_and(|p| (id, start) < p) {
                    bail!("BCF output must be coordinate-sorted in header dictionary order");
                }
                let chunk_start = writer.get_ref().virtual_position();
                writer.write_variant_record(header, record)?;
                let chunk = Chunk::new(chunk_start, writer.get_ref().virtual_position());
                indexer.add_record(Some((id, start, end, true)), chunk)?;
                *previous = Some((id, start));
            }
        }
        Ok(())
    }

    fn finish(self, path: &Path) -> Result<()> {
        match self {
            Self::Vcf(mut writer, index) => {
                writer.try_finish()?;
                index.write(path)?;
            }
            Self::Bcf {
                mut writer, indexer, reference_count, ..
            } => {
                writer.try_finish()?;
                let index_file = File::options()
                    .write(true)
                    .create_new(true)
                    .open(csi_index_path(path))?;
                let mut index_writer = csi::io::Writer::new(index_file);
                index_writer.write_index(&indexer.build(reference_count))?;
                index_writer.get_mut().try_finish()?;
            }
        }
        Ok(())
    }
}
