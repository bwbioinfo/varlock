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

/// A compact sequence dictionary entry (name + length) extracted from BAM or FASTA index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SeqEntry {
    pub(crate) name: String,
    pub(crate) length: u64,
}

/// Require identical contig names and lengths, independent of FASTA order.
/// FAI has no M5 checksums, so this checks SN/LN only, not sequence identity.
pub(crate) fn validate_bam_fasta_dict(
    bam_seqs: &[SeqEntry],
    fasta_seqs: &HashMap<String, u64>,
    fasta_path: &Path,
) -> Result<()> {
    let mut missing = Vec::new();
    let mut mismatched = Vec::new();
    let bam_names: std::collections::HashSet<_> =
        bam_seqs.iter().map(|entry| entry.name.as_str()).collect();
    let mut extra: Vec<_> = fasta_seqs
        .keys()
        .filter(|name| !bam_names.contains(name.as_str()))
        .cloned()
        .collect();
    extra.sort();

    for entry in bam_seqs {
        match fasta_seqs.get(&entry.name) {
            None => missing.push(entry.name.clone()),
            Some(&fasta_len) if fasta_len != entry.length => {
                mismatched.push(format!(
                    "{} (BAM={} FASTA={})",
                    entry.name, entry.length, fasta_len
                ));
            }
            _ => {}
        }
    }

    if missing.is_empty() && extra.is_empty() && mismatched.is_empty() {
        return Ok(());
    }

    let mut msg = format!(
        "BAM/FASTA sequence dictionary mismatch (reference: {}):",
        fasta_path.display()
    );
    for (label, entries) in [
        ("missing contigs in FASTA", missing),
        ("extra contigs in FASTA", extra),
        ("length mismatches", mismatched),
    ] {
        if !entries.is_empty() {
            let preview = entries.iter().take(5).cloned().collect::<Vec<_>>().join(", ");
            msg.push_str(&format!("\n  {label} ({}): {preview}", entries.len()));
            if entries.len() > 5 {
                msg.push_str(&format!(", ... (+{} more)", entries.len() - 5));
            }
        }
    }
    bail!("{msg}\nUse the BAM alignment reference build with matching contigs and lengths, and regenerate its .fai if stale.");
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
    /// Return a map of contig name → length derived from the FASTA index.
    pub(crate) fn fasta_lengths(&self) -> HashMap<String, u64> {
        self.records
            .iter()
            .map(|(name, rec)| (name.clone(), rec.length))
            .collect()
    }

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

    pub(crate) fn fetch_bases(
        &mut self,
        ref_name: &str,
        start_pos1: u32,
        len: u32,
    ) -> Result<Vec<u8>> {
        let mut bases = Vec::with_capacity(len as usize);
        for offset in 0..len {
            let position = start_pos1
                .checked_add(offset)
                .context("reference position overflow while fetching allele")?;
            bases.push(self.fetch_base(ref_name, position)?);
        }
        Ok(bases)
    }
}
