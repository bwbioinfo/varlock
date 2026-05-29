use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;

#[derive(Clone, Debug)]
struct FaiRecord {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

pub(crate) struct FastaIndex {
    reference_path: PathBuf,
    records: HashMap<String, FaiRecord>,
    reader: FastaReader,
    line_cache: Option<FastaLineCache>,
}

enum FastaReader {
    Plain(File),
    Bgzf(bgzf::io::IndexedReader<File>),
}

struct FastaLineCache {
    ref_name: String,
    line: u64,
    bases: Vec<u8>,
}

pub(crate) fn open_fasta_index(reference: &Path) -> Result<FastaIndex> {
    let records = open_fai(reference)?;
    let reader = if is_bgzf_reference(reference) {
        let indexed_reader = bgzf::io::indexed_reader::Builder::default()
            .build_from_path(reference)
            .with_context(|| {
                format!(
                    "failed to open indexed BGZF reference {} (expected .gzi alongside file)",
                    reference.display()
                )
            })?;
        FastaReader::Bgzf(indexed_reader)
    } else {
        let file = File::open(reference)
            .with_context(|| format!("failed to open reference {}", reference.display()))?;
        FastaReader::Plain(file)
    };

    Ok(FastaIndex {
        reference_path: reference.to_path_buf(),
        records,
        reader,
        line_cache: None,
    })
}

fn open_fai(reference: &Path) -> Result<HashMap<String, FaiRecord>> {
    let mut fai_path = reference.to_path_buf();
    if let Some(ext) = reference.extension().and_then(|ext| ext.to_str()) {
        fai_path.set_extension(format!("{}.fai", ext));
    } else {
        fai_path.set_extension("fai");
    }
    if !fai_path.exists() {
        let alt = reference.with_extension("fai");
        if alt.exists() {
            fai_path = alt;
        }
    }
    let file = File::open(&fai_path)
        .with_context(|| format!("failed to open FASTA index {}", fai_path.display()))?;
    let reader = BufReader::new(file);
    let mut records = HashMap::new();
    for line in reader.lines() {
        let line = line.context("failed to read FASTA index line")?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        let record = FaiRecord {
            name: parts[0].to_string(),
            length: parts[1].parse::<u64>().context("invalid FAI length")?,
            offset: parts[2].parse::<u64>().context("invalid FAI offset")?,
            line_bases: parts[3].parse::<u64>().context("invalid FAI line_bases")?,
            line_width: parts[4].parse::<u64>().context("invalid FAI line_width")?,
        };
        records.insert(record.name.clone(), record);
    }
    Ok(records)
}

fn is_bgzf_reference(reference: &Path) -> bool {
    reference
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "gz" | "bgz" | "bgzf"))
        .unwrap_or(false)
}

impl FastaIndex {
    pub(crate) fn reference_lengths(&self, ref_names: &[String]) -> Result<Vec<u64>> {
        let mut lengths = Vec::with_capacity(ref_names.len());
        for name in ref_names {
            let record = self
                .records
                .get(name)
                .with_context(|| format!("reference {} missing in FASTA index", name))?;
            lengths.push(record.length);
        }
        Ok(lengths)
    }

    pub(crate) fn fetch_base(&mut self, ref_name: &str, pos1: u32) -> Result<u8> {
        let record = self
            .records
            .get(ref_name)
            .with_context(|| format!("reference {} missing in FASTA index", ref_name))?;
        let length = record.length;
        let offset = record.offset;
        let line_bases = record.line_bases;
        let line_width = record.line_width;
        if pos1 == 0 || pos1 as u64 > length {
            bail!("reference position out of bounds: {}:{}", ref_name, pos1);
        }
        if line_bases == 0 {
            bail!(
                "invalid FASTA index entry for {}: line_bases cannot be 0",
                ref_name
            );
        }

        let pos0 = pos1 as u64 - 1;
        let line = pos0 / line_bases;
        let column = pos0 % line_bases;
        let column_usize = usize::try_from(column).context("FASTA line column overflow")?;

        if let Some(cache) = &self.line_cache
            && cache.ref_name == ref_name
            && cache.line == line
            && let Some(base) = cache.bases.get(column_usize).copied()
        {
            return Ok(base.to_ascii_uppercase());
        }

        let line_start_offset = offset + line * line_width;
        let bases_in_line_u64 = (length - line * line_bases).min(line_bases);
        let bases_in_line =
            usize::try_from(bases_in_line_u64).context("FASTA line length overflow")?;
        let mut line_buf = vec![0u8; bases_in_line];

        match &mut self.reader {
            FastaReader::Plain(file) => {
                file.seek(SeekFrom::Start(line_start_offset))
                    .with_context(|| {
                        format!("failed to seek reference {}", self.reference_path.display())
                    })?;
                file.read_exact(&mut line_buf).with_context(|| {
                    format!("failed to read reference {}", self.reference_path.display())
                })?
            }
            FastaReader::Bgzf(reader) => {
                reader
                    .seek(SeekFrom::Start(line_start_offset))
                    .with_context(|| {
                        format!(
                            "failed to seek indexed BGZF reference {}",
                            self.reference_path.display()
                        )
                    })?;
                reader.read_exact(&mut line_buf).with_context(|| {
                    format!(
                        "failed to read indexed BGZF reference {}",
                        self.reference_path.display()
                    )
                })?
            }
        };

        let base = *line_buf.get(column_usize).with_context(|| {
            format!(
                "reference position out of cached line bounds: {}:{}",
                ref_name, pos1
            )
        })?;
        self.line_cache = Some(FastaLineCache {
            ref_name: ref_name.to_string(),
            line,
            bases: line_buf,
        });
        Ok(base.to_ascii_uppercase())
    }
}
