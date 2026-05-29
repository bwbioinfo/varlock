use std::{
    collections::HashSet,
    ffi::{OsStr, OsString},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
    thread,
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, Result, anyhow, bail};
use flate2::read::MultiGzDecoder;

const WRAP_WIDTH: usize = 60;

#[derive(Clone, Debug)]
struct SpoolRecord {
    name: String,
    temp_path: PathBuf,
    length: u64,
}

struct ActiveRecord {
    name: String,
    temp_path: PathBuf,
    length: u64,
    writer: BufWriter<File>,
}

struct TempWorkDir {
    path: PathBuf,
}

impl Drop for TempWorkDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

pub fn prepare_reference(reference: &Path, verbose: u8, threads: usize) -> Result<PathBuf> {
    let compressed = is_compressed(reference);
    let has_fai = existing_fai_path(reference).is_some();
    let has_gzi = existing_gzi_path(reference).is_some();

    if has_fai && (!compressed || has_gzi) {
        if verbose > 1 {
            eprintln!(
                "[fasta_prep] using existing indexed reference {}",
                reference.display()
            );
        }
        return Ok(reference.to_path_buf());
    }

    let prepared = prepared_fasta_path(reference);
    if prepared.exists()
        && existing_fai_path(&prepared).is_some()
        && existing_gzi_path(&prepared).is_some()
    {
        if verbose > 0 {
            eprintln!(
                "[fasta_prep] using existing prepared reference {}",
                prepared.display()
            );
        }
        return Ok(prepared);
    }

    if verbose > 0 {
        eprintln!(
            "[fasta_prep] preparing sorted/indexed reference from {}",
            reference.display()
        );
    }
    let record_count = sort_and_index_fasta(reference, &prepared, threads, verbose)?;
    if verbose > 0 {
        eprintln!(
            "[fasta_prep] prepared records={} output={} index={} gzi={}",
            record_count,
            prepared.display(),
            primary_fai_path(&prepared).display(),
            gzi_path(&prepared).display()
        );
    }

    Ok(prepared)
}

fn sort_and_index_fasta(input: &Path, output: &Path, threads: usize, verbose: u8) -> Result<usize> {
    let work_dir = create_temp_work_dir()?;
    let mut records = parse_and_spool_records(input, &work_dir.path)?;
    records.sort_by(|a, b| a.name.cmp(&b.name));
    records = wrap_records(records, &work_dir.path, threads.max(1), verbose)?;

    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let output_file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let worker_count = NonZeroUsize::new(threads.max(1)).expect("non-zero");
    let mut writer =
        noodles_bgzf::io::MultithreadedWriter::with_worker_count(worker_count, output_file);

    let fai_path = primary_fai_path(output);
    let fai_file = File::create(&fai_path)
        .with_context(|| format!("failed to create {}", fai_path.display()))?;
    let mut fai_writer = BufWriter::new(fai_file);

    let mut offset = 0u64;
    for record in &records {
        let header_line = format!(">{}\n", record.name);
        writer
            .write_all(header_line.as_bytes())
            .with_context(|| format!("failed to write FASTA header for {}", record.name))?;
        offset += header_line.len() as u64;
        let sequence_offset = offset;

        let copied = copy_wrapped_sequence(&record.temp_path, &mut writer)
            .with_context(|| format!("failed to write FASTA sequence for {}", record.name))?;
        offset += copied;

        let (line_bases, line_width) = line_layout(record.length);
        writeln!(
            fai_writer,
            "{}\t{}\t{}\t{}\t{}",
            record.name, record.length, sequence_offset, line_bases, line_width
        )
        .with_context(|| format!("failed to write FAI line for {}", record.name))?;
    }

    writer
        .finish()
        .with_context(|| format!("failed to finish BGZF {}", output.display()))?;
    fai_writer
        .flush()
        .with_context(|| format!("failed to flush {}", fai_path.display()))?;
    build_gzi_index_for_bgzf(output)?;

    Ok(records.len())
}

fn parse_and_spool_records(input: &Path, spool_dir: &Path) -> Result<Vec<SpoolRecord>> {
    let input_file =
        File::open(input).with_context(|| format!("failed to open {}", input.display()))?;
    let mut reader: Box<dyn BufRead> = if is_compressed(input) {
        Box::new(BufReader::new(MultiGzDecoder::new(input_file)))
    } else {
        Box::new(BufReader::new(input_file))
    };

    let mut records = Vec::new();
    let mut active: Option<ActiveRecord> = None;
    let mut seen_names = HashSet::new();

    for (line_idx, line_result) in reader.by_ref().lines().enumerate() {
        let line_number = line_idx + 1;
        let line =
            line_result.with_context(|| format!("failed to read FASTA line {}", line_number))?;

        if let Some(header) = line.strip_prefix('>') {
            if let Some(current) = active.take() {
                records.push(finish_record(current)?);
            }

            let header = header.trim();
            if header.is_empty() {
                bail!("empty FASTA header at line {}", line_number);
            }
            let name = header
                .split_whitespace()
                .next()
                .context("missing FASTA sequence name")?
                .to_string();

            if !seen_names.insert(name.clone()) {
                bail!("duplicate FASTA sequence name {}", name);
            }

            let temp_path = spool_dir.join(format!("record_{:06}.seq", records.len()));
            let temp_file = File::create(&temp_path)
                .with_context(|| format!("failed to create {}", temp_path.display()))?;
            active = Some(ActiveRecord {
                name,
                temp_path,
                length: 0,
                writer: BufWriter::new(temp_file),
            });
            continue;
        }

        let cleaned = line
            .as_bytes()
            .iter()
            .copied()
            .filter(|b| !b.is_ascii_whitespace())
            .collect::<Vec<_>>();
        if cleaned.is_empty() {
            continue;
        }

        let current = active
            .as_mut()
            .with_context(|| format!("FASTA sequence before header at line {}", line_number))?;
        current
            .writer
            .write_all(&cleaned)
            .with_context(|| format!("failed to write temporary sequence for {}", current.name))?;
        current.length += cleaned.len() as u64;
    }

    if let Some(current) = active.take() {
        records.push(finish_record(current)?);
    }

    if records.is_empty() {
        bail!("no FASTA records found in {}", input.display());
    }

    Ok(records)
}

fn finish_record(mut active: ActiveRecord) -> Result<SpoolRecord> {
    active
        .writer
        .flush()
        .with_context(|| format!("failed to flush temporary sequence for {}", active.name))?;
    Ok(SpoolRecord {
        name: active.name,
        temp_path: active.temp_path,
        length: active.length,
    })
}

fn wrap_records(
    mut records: Vec<SpoolRecord>,
    work_dir: &Path,
    threads: usize,
    verbose: u8,
) -> Result<Vec<SpoolRecord>> {
    let worker_count = threads.min(records.len().max(1));
    if verbose > 0 {
        eprintln!(
            "[fasta_prep] wrapping {} record(s) with {} thread(s)",
            records.len(),
            worker_count
        );
    }

    if records.is_empty() {
        return Ok(records);
    }

    if worker_count <= 1 {
        for (idx, record) in records.iter_mut().enumerate() {
            let wrapped = work_dir.join(format!("wrapped_{idx:06}.seq"));
            wrap_sequence_file(&record.temp_path, &wrapped)?;
            let _ = std::fs::remove_file(&record.temp_path);
            record.temp_path = wrapped;
        }
        return Ok(records);
    }

    let mut handles = Vec::with_capacity(worker_count);
    for worker_id in 0..worker_count {
        let assigned = records
            .iter()
            .cloned()
            .enumerate()
            .filter(|(idx, _)| idx % worker_count == worker_id)
            .collect::<Vec<_>>();
        let work_dir = work_dir.to_path_buf();

        handles.push(thread::spawn(move || -> Result<Vec<(usize, PathBuf)>> {
            let mut out = Vec::with_capacity(assigned.len());
            for (idx, record) in assigned {
                let wrapped = work_dir.join(format!("wrapped_{idx:06}.seq"));
                wrap_sequence_file(&record.temp_path, &wrapped)?;
                let _ = std::fs::remove_file(&record.temp_path);
                out.push((idx, wrapped));
            }
            Ok(out)
        }));
    }

    let mut wrapped_paths = std::iter::repeat_with(|| None)
        .take(records.len())
        .collect::<Vec<Option<PathBuf>>>();
    for handle in handles {
        let worker_results = handle
            .join()
            .map_err(|_| anyhow!("fasta_prep worker thread panicked"))??;
        for (idx, wrapped) in worker_results {
            wrapped_paths[idx] = Some(wrapped);
        }
    }

    for (idx, record) in records.iter_mut().enumerate() {
        record.temp_path = wrapped_paths[idx]
            .take()
            .with_context(|| format!("missing wrapped sequence for record {}", record.name))?;
    }

    Ok(records)
}

fn wrap_sequence_file(input: &Path, output: &Path) -> Result<()> {
    let input_file =
        File::open(input).with_context(|| format!("failed to open {}", input.display()))?;
    let mut reader = BufReader::new(input_file);
    let output_file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::new(output_file);

    let mut buf = [0u8; 8192];
    let mut col = 0usize;
    let mut wrote_any = false;

    loop {
        let n = reader
            .read(&mut buf)
            .with_context(|| format!("failed to read {}", input.display()))?;
        if n == 0 {
            break;
        }
        wrote_any = true;
        for &b in &buf[..n] {
            writer
                .write_all(&[b])
                .with_context(|| format!("failed to write {}", output.display()))?;
            col += 1;
            if col == WRAP_WIDTH {
                writer
                    .write_all(b"\n")
                    .with_context(|| format!("failed to write {}", output.display()))?;
                col = 0;
            }
        }
    }

    if !wrote_any || col != 0 {
        writer
            .write_all(b"\n")
            .with_context(|| format!("failed to write {}", output.display()))?;
    }

    writer
        .flush()
        .with_context(|| format!("failed to flush {}", output.display()))?;
    Ok(())
}

fn copy_wrapped_sequence<W: Write>(sequence_path: &Path, writer: &mut W) -> Result<u64> {
    let file = File::open(sequence_path)
        .with_context(|| format!("failed to open {}", sequence_path.display()))?;
    let mut reader = BufReader::new(file);
    std::io::copy(&mut reader, writer)
        .with_context(|| format!("failed to copy {}", sequence_path.display()))
}

fn line_layout(length: u64) -> (u64, u64) {
    if length == 0 {
        (1, 2)
    } else {
        let line_bases = length.min(WRAP_WIDTH as u64);
        (line_bases, line_bases + 1)
    }
}

fn create_temp_work_dir() -> Result<TempWorkDir> {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system clock before UNIX_EPOCH")?
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "varlock_fasta_prep_{}_{}",
        std::process::id(),
        now
    ));
    std::fs::create_dir_all(&path)
        .with_context(|| format!("failed to create temporary directory {}", path.display()))?;
    Ok(TempWorkDir { path })
}

fn prepared_fasta_path(reference: &Path) -> PathBuf {
    let parent = reference.parent().unwrap_or_else(|| Path::new("."));
    let filename = reference
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or("reference.fa");
    let base = filename
        .strip_suffix(".fa.gz")
        .or_else(|| filename.strip_suffix(".fasta.gz"))
        .or_else(|| filename.strip_suffix(".fna.gz"))
        .or_else(|| filename.strip_suffix(".gz"))
        .or_else(|| filename.strip_suffix(".fa"))
        .or_else(|| filename.strip_suffix(".fasta"))
        .or_else(|| filename.strip_suffix(".fna"))
        .unwrap_or(filename);
    parent.join(format!("{base}.sorted.fa.gz"))
}

fn is_compressed(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "gz" | "bgz" | "bgzf"))
        .unwrap_or(false)
}

fn existing_fai_path(reference: &Path) -> Option<PathBuf> {
    let primary = primary_fai_path(reference);
    if primary.exists() {
        return Some(primary);
    }
    let alternate = reference.with_extension("fai");
    if alternate.exists() {
        return Some(alternate);
    }
    None
}

fn existing_gzi_path(reference: &Path) -> Option<PathBuf> {
    let path = gzi_path(reference);
    if path.exists() { Some(path) } else { None }
}

fn primary_fai_path(reference: &Path) -> PathBuf {
    let mut fai_path = reference.to_path_buf();
    if let Some(ext) = reference.extension().and_then(|ext| ext.to_str()) {
        fai_path.set_extension(format!("{}.fai", ext));
    } else {
        fai_path.set_extension("fai");
    }
    fai_path
}

fn gzi_path(reference: &Path) -> PathBuf {
    push_ext(reference.to_path_buf(), "gzi")
}

fn push_ext(path: PathBuf, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut s = OsString::from(path);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn build_gzi_index_for_bgzf(path: &Path) -> Result<()> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = noodles_bgzf::io::Reader::new(file);

    let mut entries = Vec::new();
    let mut prev_compressed = reader.virtual_position().compressed();
    let mut uncompressed_pos = 0u64;

    loop {
        let compressed = reader.virtual_position().compressed();
        let buf = reader
            .fill_buf()
            .with_context(|| format!("failed to scan BGZF {}", path.display()))?;
        if buf.is_empty() {
            break;
        }

        if compressed != prev_compressed {
            entries.push((compressed, uncompressed_pos));
            prev_compressed = compressed;
        }

        let len = buf.len();
        reader.consume(len);
        uncompressed_pos += len as u64;
    }

    let index = noodles_bgzf::gzi::Index::from(entries);
    let gzi = gzi_path(path);
    noodles_bgzf::gzi::fs::write(&gzi, &index)
        .with_context(|| format!("failed to write {}", gzi.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Read, io::Write, path::Path};
    use anyhow::{Context, Result};
    use flate2::{Compression, write::GzEncoder};
    use tempfile::tempdir;
    use super::{gzi_path, prepare_reference, primary_fai_path};

    fn read_bgzf_to_string(path: &Path) -> Result<String> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        let mut reader = noodles_bgzf::io::Reader::new(file);
        let mut s = String::new();
        reader
            .read_to_string(&mut s)
            .with_context(|| format!("failed to read {}", path.display()))?;
        Ok(s)
    }

    #[test]
    fn prepares_sorted_and_indexed_reference_when_unindexed() -> Result<()> {
        let dir = tempdir().context("failed to create tempdir")?;
        let reference = dir.path().join("ref.fa");
        std::fs::write(&reference, b">chr2\nTTAA\n>chr1\nACGT\n")?;

        let prepared = prepare_reference(&reference, 0, 4)?;
        assert_ne!(prepared, reference);
        assert!(prepared.exists());
        assert!(primary_fai_path(&prepared).exists());
        assert!(gzi_path(&prepared).exists());
        assert!(prepared
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("gz"))
            .unwrap_or(false));

        let text = read_bgzf_to_string(&prepared)?;
        assert!(text.starts_with(">chr1\n"), "expected chr1 first, got: {text}");
        Ok(())
    }

    #[test]
    fn prepares_from_gz_reference() -> Result<()> {
        let dir = tempdir()?;
        let reference = dir.path().join("ref.fa.gz");
        let file = File::create(&reference)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(b">chrB\nTTAA\n>chrA\nACGT\n")?;
        encoder.finish()?;

        let prepared = prepare_reference(&reference, 0, 4)?;
        assert!(prepared.exists());
        assert!(primary_fai_path(&prepared).exists());
        assert!(gzi_path(&prepared).exists());

        let text = read_bgzf_to_string(&prepared)?;
        assert!(text.starts_with(">chrA\n"), "expected chrA first, got: {text}");
        Ok(())
    }

    #[test]
    fn returns_original_when_reference_already_indexed() -> Result<()> {
        let dir = tempdir()?;
        let reference = dir.path().join("ref.fa");
        std::fs::write(&reference, b">chr1\nACGT\n")?;
        std::fs::write(primary_fai_path(&reference), b"chr1\t4\t6\t4\t5\n")?;

        let prepared = prepare_reference(&reference, 0, 2)?;
        assert_eq!(prepared, reference);
        Ok(())
    }

    #[test]
    fn reuses_existing_prepared_reference() -> Result<()> {
        let dir = tempdir()?;
        let reference = dir.path().join("ref.fa");
        std::fs::write(&reference, b">chr1\nACGT\n")?;

        let first = prepare_reference(&reference, 0, 2)?;
        assert_ne!(first, reference);

        // Second call should return the same prepared path without reprocessing
        let second = prepare_reference(&reference, 0, 2)?;
        assert_eq!(first, second);
        Ok(())
    }
}
