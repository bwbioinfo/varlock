use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::{
    self,
    index::reference_sequence::bin::Chunk,
    index::reference_sequence::index::BinnedIndex,
    index::{
        Header as TabixHeader,
        header::{Format as TabixFormat, ReferenceSequenceNames},
    },
};
use noodles_tabix as tabix;

use crate::{ExecutionContext, FilterArgs, IndexType, log_verbose};

#[derive(Clone, Debug)]
struct MaxInfoFilter {
    field: String,
    max: f64,
}

#[derive(Debug)]
struct FilterSpec {
    require_info: HashSet<String>,
    exclude_info: HashSet<String>,
    max_info: Vec<MaxInfoFilter>,
}

pub(crate) fn run(args: FilterArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    let spec = FilterSpec::from_args(&args)?;
    let metrics = filter_vcf(&args.input, &args.output, args.index_type, &spec)?;

    log_verbose(
        ctx,
        format!(
            "filter stage=done input_records={} output_records={} output={} elapsed={:.2?}",
            metrics.input_records,
            metrics.output_records,
            args.output.display(),
            started.elapsed()
        ),
    );
    Ok(())
}

impl FilterSpec {
    fn from_args(args: &FilterArgs) -> Result<Self> {
        let max_info = args
            .max_info
            .iter()
            .map(|value| {
                let (field, max) = value.split_once('=').with_context(|| {
                    format!("invalid --max-info {value:?}; expected FIELD=VALUE")
                })?;
                if field.is_empty() || max.is_empty() {
                    bail!("invalid --max-info {value:?}; expected FIELD=VALUE");
                }
                let max = max
                    .parse::<f64>()
                    .with_context(|| format!("invalid --max-info threshold {max:?}"))?;
                Ok(MaxInfoFilter {
                    field: field.to_string(),
                    max,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            require_info: args.require_info.iter().cloned().collect(),
            exclude_info: args.exclude_info.iter().cloned().collect(),
            max_info,
        })
    }

    fn keep_record(&self, info_text: &str) -> Result<bool> {
        let info = parse_info(info_text);
        for field in &self.require_info {
            if !info.contains_key(field.as_str()) {
                return Ok(false);
            }
        }
        for field in &self.exclude_info {
            if info.contains_key(field.as_str()) {
                return Ok(false);
            }
        }
        for filter in &self.max_info {
            let Some(value) = info.get(filter.field.as_str()) else {
                return Ok(false);
            };
            if !all_numeric_values_at_most(value, filter.max)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

#[derive(Default)]
struct FilterMetrics {
    input_records: usize,
    output_records: usize,
}

fn filter_vcf(
    input: &Path,
    output: &Path,
    index_type: IndexType,
    spec: &FilterSpec,
) -> Result<FilterMetrics> {
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

    let mut line = String::new();
    let mut saw_column_header = false;
    let mut metrics = FilterMetrics::default();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            saw_column_header = true;
            writeln!(writer, "{trimmed}")?;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.input_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            if spec.keep_record(fields[7])? {
                let record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{trimmed}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }
    if !saw_column_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

#[derive(Clone, Debug)]
struct IndexRecord {
    chrom: String,
    position: Position,
}

impl IndexRecord {
    fn from_fields(fields: &[&str]) -> Result<Self> {
        let pos = fields[1]
            .parse::<usize>()
            .with_context(|| format!("invalid VCF position {:?}", fields[1]))?;
        let position = Position::try_from(pos).context("invalid VCF position for indexing")?;
        Ok(Self {
            chrom: fields[0].to_string(),
            position,
        })
    }
}

enum OutputIndex {
    Csi {
        indexer: binning_index::Indexer<BinnedIndex>,
        reference_ids: HashMap<String, usize>,
        reference_names: Vec<String>,
    },
    Tbi(tabix::index::Indexer),
}

impl OutputIndex {
    fn new(index_type: IndexType) -> Self {
        match index_type {
            IndexType::Csi => Self::Csi {
                indexer: binning_index::Indexer::<BinnedIndex>::default(),
                reference_ids: HashMap::new(),
                reference_names: Vec::new(),
            },
            IndexType::Tbi => {
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
                Self::Tbi(indexer)
            }
        }
    }

    fn add_record(&mut self, record: &IndexRecord, chunk: Chunk) -> Result<()> {
        match self {
            Self::Csi {
                indexer,
                reference_ids,
                reference_names,
            } => {
                let reference_sequence_id =
                    reference_id_for(reference_ids, reference_names, &record.chrom);
                indexer
                    .add_record(
                        Some((
                            reference_sequence_id,
                            record.position,
                            record.position,
                            true,
                        )),
                        chunk,
                    )
                    .context("failed to update CSI index")?;
            }
            Self::Tbi(indexer) => indexer
                .add_record(&record.chrom, record.position, record.position, chunk)
                .context("failed to update TBI index")?,
        }
        Ok(())
    }

    fn write(self, output: &Path) -> Result<()> {
        match self {
            Self::Csi {
                indexer,
                reference_names,
                ..
            } => {
                let mut csi_reference_names = ReferenceSequenceNames::new();
                for name in &reference_names {
                    csi_reference_names.insert(name.as_str().into());
                }
                let header = TabixHeader::builder()
                    .set_format(TabixFormat::Vcf)
                    .set_reference_sequence_name_index(0)
                    .set_start_position_index(1)
                    .set_end_position_index(None)
                    .set_line_comment_prefix(b'#')
                    .set_line_skip_count(0)
                    .set_reference_sequence_names(csi_reference_names)
                    .build();
                let index = indexer.set_header(header).build(reference_names.len());
                let index_path = index_path(output, "csi")?;
                let index_file = File::create(&index_path)
                    .with_context(|| format!("failed to create index {}", index_path.display()))?;
                let mut writer = noodles_csi::io::Writer::new(index_file);
                writer
                    .write_index(&index)
                    .context("failed to write CSI index")?;
            }
            Self::Tbi(indexer) => {
                let index = indexer.build();
                let index_path = index_path(output, "tbi")?;
                let index_file = File::create(&index_path)
                    .with_context(|| format!("failed to create index {}", index_path.display()))?;
                let mut writer = tabix::io::Writer::new(index_file);
                writer
                    .write_index(&index)
                    .context("failed to write TBI index")?;
            }
        }
        Ok(())
    }
}

fn reference_id_for(
    reference_ids: &mut HashMap<String, usize>,
    reference_names: &mut Vec<String>,
    chrom: &str,
) -> usize {
    if let Some(id) = reference_ids.get(chrom) {
        *id
    } else {
        let id = reference_names.len();
        reference_ids.insert(chrom.to_string(), id);
        reference_names.push(chrom.to_string());
        id
    }
}

fn parse_info(info: &str) -> HashMap<&str, &str> {
    let mut out = HashMap::new();
    if info == "." || info.is_empty() {
        return out;
    }
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            out.insert(key, value);
        } else {
            out.insert(item, "1");
        }
    }
    out
}

fn all_numeric_values_at_most(value: &str, max: f64) -> Result<bool> {
    let mut saw_value = false;
    for part in value.split(',') {
        if part == "." || part.is_empty() {
            continue;
        }
        saw_value = true;
        let parsed = part
            .parse::<f64>()
            .with_context(|| format!("INFO value {part:?} is not numeric"))?;
        if parsed > max {
            return Ok(false);
        }
    }
    Ok(saw_value)
}

fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if is_gz_path(path) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn is_gz_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "gz" | "bgz" | "bgzf"))
        .unwrap_or(false)
}

fn index_path(output: &Path, ext: &str) -> Result<PathBuf> {
    output
        .file_name()
        .and_then(|name| name.to_str())
        .map(|name| output.with_file_name(format!("{name}.{ext}")))
        .context("invalid output path")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn filter_spec_requires_and_excludes_info() -> Result<()> {
        let spec = FilterSpec {
            require_info: ["AF".to_string()].into_iter().collect(),
            exclude_info: ["COMMON".to_string()].into_iter().collect(),
            max_info: Vec::new(),
        };
        assert!(spec.keep_record("AF=0.1")?);
        assert!(!spec.keep_record("DP=10")?);
        assert!(!spec.keep_record("AF=0.1;COMMON")?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_max_info_to_all_numeric_values() -> Result<()> {
        let spec = FilterSpec {
            require_info: HashSet::new(),
            exclude_info: HashSet::new(),
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
        };
        assert!(spec.keep_record("AF=0.01,0.1")?);
        assert!(!spec.keep_record("AF=0.01,0.2")?);
        assert!(!spec.keep_record("DP=10")?);
        Ok(())
    }

    #[test]
    fn filter_vcf_writes_filtered_bgzip_and_index() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let output = dir.path().join("out.vcf.gz");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.2\n",
        )?;
        let spec = FilterSpec {
            require_info: HashSet::new(),
            exclude_info: HashSet::new(),
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
        };

        let metrics = filter_vcf(&input, &output, IndexType::Csi, &spec)?;
        assert_eq!(metrics.input_records, 2);
        assert_eq!(metrics.output_records, 1);
        assert!(dir.path().join("out.vcf.gz.csi").exists());

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        Ok(())
    }
}
