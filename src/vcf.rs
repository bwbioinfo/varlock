use std::{
    collections::HashMap,
    fs::File,
    hash::{Hash, Hasher},
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
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

use crate::IndexType;

#[derive(Clone, Debug, Eq)]
pub(crate) struct VariantKey {
    pub(crate) chrom: String,
    pub(crate) pos: String,
    pub(crate) ref_allele: String,
    pub(crate) alt_allele: String,
}

impl PartialEq for VariantKey {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom
            && self.pos == other.pos
            && self.ref_allele == other.ref_allele
            && self.alt_allele == other.alt_allele
    }
}

impl Hash for VariantKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chrom.hash(state);
        self.pos.hash(state);
        self.ref_allele.hash(state);
        self.alt_allele.hash(state);
    }
}

#[derive(Clone, Debug)]
pub(crate) struct IndexRecord {
    pub(crate) chrom: String,
    pub(crate) position: Position,
}

impl IndexRecord {
    pub(crate) fn from_fields(fields: &[&str]) -> Result<Self> {
        let pos = fields
            .get(1)
            .context("invalid VCF record missing POS")?
            .parse::<usize>()
            .with_context(|| format!("invalid VCF position {:?}", fields[1]))?;
        let position = Position::try_from(pos).context("invalid VCF position for indexing")?;
        let chrom = fields
            .first()
            .context("invalid VCF record missing CHROM")?
            .to_string();
        Ok(Self { chrom, position })
    }

    pub(crate) fn from_line(line: &str) -> Result<Self> {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 8 {
            bail!("invalid VCF record with fewer than 8 fields: {line}");
        }
        Self::from_fields(&fields)
    }
}

pub(crate) struct VcfRecord<'a> {
    pub(crate) info_text: &'a str,
    fields: &'a [&'a str],
    sample_index: &'a HashMap<String, usize>,
    format_index: HashMap<&'a str, usize>,
}

impl<'a> VcfRecord<'a> {
    pub(crate) fn new(fields: &'a [&'a str], sample_index: &'a HashMap<String, usize>) -> Self {
        let format_index = if fields.len() >= 9 {
            fields[8]
                .split(':')
                .enumerate()
                .map(|(i, field)| (field, i))
                .collect()
        } else {
            HashMap::new()
        };
        Self {
            info_text: fields[7],
            fields,
            sample_index,
            format_index,
        }
    }

    pub(crate) fn sample_field(&self, sample: &str, field: &str) -> Option<&'a str> {
        let sample_offset = *self.sample_index.get(sample)?;
        let format_offset = *self.format_index.get(field)?;
        let sample_text = self.fields.get(9 + sample_offset)?;
        sample_text.split(':').nth(format_offset)
    }
}

pub(crate) fn parse_header_samples(header_line: &str) -> (Vec<String>, HashMap<String, usize>) {
    let fields = header_line.split('\t').collect::<Vec<_>>();
    if fields.len() <= 9 {
        return (Vec::new(), HashMap::new());
    }
    let sample_names = fields[9..]
        .iter()
        .map(|sample| (*sample).to_string())
        .collect::<Vec<_>>();
    let sample_index = sample_names
        .iter()
        .enumerate()
        .map(|(i, sample)| (sample.clone(), i))
        .collect();
    (sample_names, sample_index)
}

pub(crate) fn variant_keys_from_fields(fields: &[&str]) -> Vec<VariantKey> {
    let Some(chrom) = fields.first() else {
        return Vec::new();
    };
    let Some(pos) = fields.get(1) else {
        return Vec::new();
    };
    let Some(ref_allele) = fields.get(3) else {
        return Vec::new();
    };
    let Some(alt_field) = fields.get(4) else {
        return Vec::new();
    };
    if *alt_field == "." || alt_field.is_empty() {
        return Vec::new();
    }

    alt_field
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .map(|alt| VariantKey {
            chrom: (*chrom).to_string(),
            pos: (*pos).to_string(),
            ref_allele: (*ref_allele).to_string(),
            alt_allele: alt.to_string(),
        })
        .collect()
}

pub(crate) fn parse_info_string(info: &str) -> HashMap<&str, String> {
    let mut out = HashMap::new();
    if info == "." || info.is_empty() {
        return out;
    }
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            out.insert(key, value.to_string());
        } else {
            out.insert(item, "1".to_string());
        }
    }
    out
}

pub(crate) fn parse_info_ref(info: &str) -> HashMap<&str, &str> {
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

pub(crate) fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
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

pub(crate) fn tabix_index_path(path: &Path) -> PathBuf {
    path.with_file_name(format!(
        "{}.tbi",
        path.file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default()
    ))
}

pub(crate) fn csi_index_path(path: &Path) -> PathBuf {
    path.with_file_name(format!(
        "{}.csi",
        path.file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default()
    ))
}

pub(crate) fn index_path(output: &Path, ext: &str) -> Result<PathBuf> {
    output
        .file_name()
        .and_then(|name| name.to_str())
        .map(|name| output.with_file_name(format!("{name}.{ext}")))
        .context("invalid output path")
}

pub(crate) enum OutputIndex {
    Csi {
        indexer: binning_index::Indexer<BinnedIndex>,
        reference_ids: HashMap<String, usize>,
        reference_names: Vec<String>,
    },
    Tbi(tabix::index::Indexer),
}

impl OutputIndex {
    pub(crate) fn new(index_type: IndexType) -> Self {
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

    pub(crate) fn add_record(&mut self, record: &IndexRecord, chunk: Chunk) -> Result<()> {
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

    pub(crate) fn write(self, output: &Path) -> Result<()> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variant_keys_split_multi_alt_records() {
        let keys = variant_keys_from_fields(&["chr1", "10", ".", "A", "C,G", ".", ".", "."]);
        assert_eq!(keys.len(), 2);
        assert_eq!(keys[0].alt_allele, "C");
        assert_eq!(keys[1].alt_allele, "G");
    }

    #[test]
    fn parse_header_samples_returns_names_and_index() {
        let (samples, index) =
            parse_header_samples("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb");
        assert_eq!(samples, vec!["a", "b"]);
        assert_eq!(index["a"], 0);
        assert_eq!(index["b"], 1);
    }

    #[test]
    fn vcf_record_reads_format_fields() {
        let sample_index = [
            ("tumor".to_string(), 0usize),
            ("normal".to_string(), 1usize),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        let fields = [
            "chr1", "10", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/1:12", "0/0:20",
        ];
        let record = VcfRecord::new(&fields, &sample_index);
        assert_eq!(record.sample_field("tumor", "GT"), Some("0/1"));
        assert_eq!(record.sample_field("normal", "DP"), Some("20"));
    }
}
