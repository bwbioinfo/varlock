use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    hash::{Hash, Hasher},
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
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

use crate::{AnnotateArgs, ExecutionContext, IndexType, log_verbose};

#[derive(Clone, Debug)]
struct DatabaseSpec {
    name: String,
    path: PathBuf,
}

#[derive(Clone, Debug)]
struct FieldMapping {
    src: String,
    dest: String,
}

#[derive(Clone, Debug, Eq)]
struct VariantKey {
    chrom: String,
    pos: String,
    ref_allele: String,
    alt_allele: String,
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

type AnnotationValues = BTreeMap<String, String>;
type AnnotationMap = HashMap<VariantKey, AnnotationValues>;

pub(crate) fn run(args: AnnotateArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    let databases = parse_database_specs(&args.databases)?;
    let mappings = parse_annotation_mappings(&args.annotations, &databases)?;
    let mut lookup = AnnotationLookup::open(&databases, &mappings, ctx)?;

    annotate_vcf(
        &args.input,
        &args.output,
        args.index_type,
        &databases,
        &mappings,
        &mut lookup,
    )?;

    log_verbose(
        ctx,
        format!(
            "annotate stage=done databases={} loaded_sites={} output={} elapsed={:.2?}",
            databases.len(),
            lookup.loaded_site_count(),
            args.output.display(),
            started.elapsed()
        ),
    );
    Ok(())
}

fn parse_database_specs(values: &[String]) -> Result<Vec<DatabaseSpec>> {
    let mut out = Vec::with_capacity(values.len());
    let mut seen = HashMap::new();
    for value in values {
        let (name, path) = value
            .split_once('=')
            .with_context(|| format!("invalid --database {value:?}; expected NAME=VCF"))?;
        if name.is_empty() || path.is_empty() {
            bail!("invalid --database {value:?}; expected NAME=VCF");
        }
        if seen.insert(name.to_string(), ()).is_some() {
            bail!("duplicate annotation database name {name:?}");
        }
        out.push(DatabaseSpec {
            name: name.to_string(),
            path: PathBuf::from(path),
        });
    }
    Ok(out)
}

fn parse_annotation_mappings(
    values: &[String],
    databases: &[DatabaseSpec],
) -> Result<HashMap<String, Vec<FieldMapping>>> {
    let known = databases
        .iter()
        .map(|db| (db.name.as_str(), ()))
        .collect::<HashMap<_, _>>();
    let mut out: HashMap<String, Vec<FieldMapping>> = HashMap::new();
    for value in values {
        let (db_name, mapping_text) = value
            .split_once(':')
            .with_context(|| format!("invalid --annotation {value:?}; expected NAME:SRC=DEST"))?;
        if !known.contains_key(db_name) {
            bail!("annotation mapping references unknown database {db_name:?}");
        }
        if mapping_text.is_empty() {
            bail!("annotation mapping for {db_name:?} is empty");
        }
        for mapping in mapping_text.split(',') {
            let (src, dest) = mapping.split_once('=').with_context(|| {
                format!("invalid annotation mapping {mapping:?}; expected SRC=DEST")
            })?;
            if src.is_empty() || dest.is_empty() {
                bail!("invalid annotation mapping {mapping:?}; expected SRC=DEST");
            }
            out.entry(db_name.to_string())
                .or_default()
                .push(FieldMapping {
                    src: src.to_string(),
                    dest: dest.to_string(),
                });
        }
    }
    Ok(out)
}

enum DatabaseLookup {
    InMemory {
        annotations: AnnotationMap,
    },
    Tabix {
        name: String,
        mappings: Vec<FieldMapping>,
        reader: noodles_csi::io::IndexedReader<bgzf::io::Reader<File>, tabix::Index>,
    },
}

struct AnnotationLookup {
    databases: Vec<DatabaseLookup>,
}

impl AnnotationLookup {
    fn open(
        databases: &[DatabaseSpec],
        mappings: &HashMap<String, Vec<FieldMapping>>,
        ctx: &ExecutionContext,
    ) -> Result<Self> {
        let mut opened = Vec::with_capacity(databases.len());
        for db in databases {
            let db_mappings = mappings.get(&db.name).with_context(|| {
                format!(
                    "no --annotation mapping provided for database {:?}",
                    db.name
                )
            })?;

            if tabix_index_path(&db.path).exists() {
                let reader = tabix::io::indexed_reader::Builder::default()
                    .build_from_path(&db.path)
                    .with_context(|| {
                        format!(
                            "failed to open tabix-indexed annotation database {}",
                            db.path.display()
                        )
                    })?;
                log_verbose(
                    ctx,
                    format!(
                        "annotate stage=open_database name={} mode=tabix path={}",
                        db.name,
                        db.path.display()
                    ),
                );
                opened.push(DatabaseLookup::Tabix {
                    name: db.name.clone(),
                    mappings: db_mappings.clone(),
                    reader,
                });
            } else {
                let annotations = load_one_annotation_database(db, db_mappings)?;
                log_verbose(
                    ctx,
                    format!(
                        "annotate stage=open_database name={} mode=in_memory path={} sites={}",
                        db.name,
                        db.path.display(),
                        annotations.len()
                    ),
                );
                opened.push(DatabaseLookup::InMemory { annotations });
            }
        }

        Ok(Self { databases: opened })
    }

    fn lookup(&mut self, key: &VariantKey) -> Result<AnnotationValues> {
        let mut out = AnnotationValues::new();
        for db in &mut self.databases {
            match db {
                DatabaseLookup::InMemory { annotations, .. } => {
                    if let Some(values) = annotations.get(key) {
                        out.extend(values.clone());
                    }
                }
                DatabaseLookup::Tabix {
                    name,
                    mappings,
                    reader,
                } => {
                    out.extend(query_tabix_database(name, mappings, reader, key)?);
                }
            }
        }
        Ok(out)
    }

    fn loaded_site_count(&self) -> usize {
        self.databases
            .iter()
            .map(|db| match db {
                DatabaseLookup::InMemory { annotations, .. } => annotations.len(),
                DatabaseLookup::Tabix { .. } => 0,
            })
            .sum()
    }
}

fn load_one_annotation_database(
    db: &DatabaseSpec,
    db_mappings: &[FieldMapping],
) -> Result<AnnotationMap> {
    let mut out: AnnotationMap = HashMap::new();
    let mut reader = open_text_reader(&db.path)
        .with_context(|| format!("failed to open annotation database {}", db.path.display()))?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && let Some((key, values)) = parse_database_record(trimmed, db_mappings)?
        {
            out.entry(key).or_default().extend(values);
        }
        line.clear();
    }
    Ok(out)
}

fn annotate_vcf(
    input: &Path,
    output: &Path,
    index_type: IndexType,
    databases: &[DatabaseSpec],
    mappings: &HashMap<String, Vec<FieldMapping>>,
    lookup: &mut AnnotationLookup,
) -> Result<()> {
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
    let mut wrote_annotation_headers = false;
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") && !wrote_annotation_headers {
            write_annotation_headers(&mut writer, databases, mappings)?;
            wrote_annotation_headers = true;
            writeln!(writer, "{trimmed}")?;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            let Some(index_record) = index_record_from_line(trimmed)? else {
                bail!("invalid VCF record for indexing: {trimmed}");
            };
            let chunk_start = writer.virtual_position();
            let annotated = annotate_record_line(trimmed, lookup)?;
            writeln!(writer, "{annotated}")?;
            let chunk_end = writer.virtual_position();
            output_index.add_record(&index_record, Chunk::new(chunk_start, chunk_end))?;
        }
        line.clear();
    }

    if !wrote_annotation_headers {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(())
}

fn write_annotation_headers<W: Write>(
    writer: &mut W,
    databases: &[DatabaseSpec],
    mappings: &HashMap<String, Vec<FieldMapping>>,
) -> Result<()> {
    for db in databases {
        writeln!(
            writer,
            "##varlock_annotation_database=<ID={},Path={}>",
            db.name,
            db.path.display()
        )?;
        if let Some(db_mappings) = mappings.get(&db.name) {
            for mapping in db_mappings {
                writeln!(
                    writer,
                    "##INFO=<ID={},Number=1,Type=String,Description=\"Annotation from {}:{}\">",
                    mapping.dest, db.name, mapping.src
                )?;
            }
        }
    }
    Ok(())
}

fn annotate_record_line(line: &str, lookup: &mut AnnotationLookup) -> Result<String> {
    let mut fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid VCF record with fewer than 8 fields: {line}");
    }
    let Some(key) = variant_key_from_fields(&fields) else {
        return Ok(line.to_string());
    };
    let values = lookup.lookup(&key)?;
    if values.is_empty() {
        return Ok(line.to_string());
    }

    let mut info = if fields[7] == "." || fields[7].is_empty() {
        String::new()
    } else {
        fields[7].to_string()
    };
    for (key, value) in values {
        if !info.is_empty() {
            info.push(';');
        }
        info.push_str(&key);
        info.push('=');
        info.push_str(&value);
    }
    fields[7] = &info;
    Ok(fields.join("\t"))
}

#[derive(Clone, Debug)]
struct IndexRecord {
    chrom: String,
    position: Position,
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

fn index_record_from_line(line: &str) -> Result<Option<IndexRecord>> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid VCF record with fewer than 8 fields: {line}");
    }
    let pos = fields[1]
        .parse::<usize>()
        .with_context(|| format!("invalid VCF position {:?}", fields[1]))?;
    let position = Position::try_from(pos).context("invalid VCF position for indexing")?;
    Ok(Some(IndexRecord {
        chrom: fields[0].to_string(),
        position,
    }))
}

fn query_tabix_database(
    db_name: &str,
    mappings: &[FieldMapping],
    reader: &mut noodles_csi::io::IndexedReader<bgzf::io::Reader<File>, tabix::Index>,
    key: &VariantKey,
) -> Result<AnnotationValues> {
    let region = format!("{}:{}-{}", key.chrom, key.pos, key.pos)
        .parse::<Region>()
        .with_context(|| {
            format!(
                "failed to build tabix query region for database {db_name:?}: {}:{}-{}",
                key.chrom, key.pos, key.pos
            )
        })?;
    let mut out = AnnotationValues::new();
    for result in reader
        .query(&region)
        .with_context(|| format!("failed to query tabix database {db_name:?}"))?
    {
        let record = result.with_context(|| format!("failed to read tabix record {db_name:?}"))?;
        let line = record.as_ref();
        if let Some((record_key, values)) = parse_database_record(line, mappings)?
            && record_key == *key
        {
            out.extend(values);
        }
    }
    Ok(out)
}

fn parse_database_record(
    line: &str,
    mappings: &[FieldMapping],
) -> Result<Option<(VariantKey, AnnotationValues)>> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid database VCF record with fewer than 8 fields: {line}");
    }
    let Some(key) = variant_key_from_fields(&fields) else {
        return Ok(None);
    };
    let info = parse_info(fields[7]);
    let mut values = AnnotationValues::new();
    for mapping in mappings {
        if let Some(value) = info.get(mapping.src.as_str()) {
            values.insert(mapping.dest.clone(), value.clone());
        }
    }
    if values.is_empty() {
        return Ok(None);
    }
    Ok(Some((key, values)))
}

fn variant_key_from_fields(fields: &[&str]) -> Option<VariantKey> {
    let alt = fields.get(4)?;
    if alt.contains(',') || *alt == "." || alt.is_empty() {
        return None;
    }
    Some(VariantKey {
        chrom: fields.first()?.to_string(),
        pos: fields.get(1)?.to_string(),
        ref_allele: fields.get(3)?.to_string(),
        alt_allele: alt.to_string(),
    })
}

fn parse_info(info: &str) -> HashMap<&str, String> {
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

fn tabix_index_path(path: &Path) -> PathBuf {
    path.with_file_name(format!(
        "{}.tbi",
        path.file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default()
    ))
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
    use anyhow::Result;
    use std::io::{Read, Write};
    use tempfile::tempdir;

    #[test]
    fn parses_database_and_annotation_specs() -> Result<()> {
        let databases = parse_database_specs(&["gnomad=db.vcf.gz".to_string()])?;
        assert_eq!(databases[0].name, "gnomad");
        assert_eq!(databases[0].path, PathBuf::from("db.vcf.gz"));

        let mappings = parse_annotation_mappings(
            &["gnomad:AF=gnomAD_AF,AC=gnomAD_AC".to_string()],
            &databases,
        )?;
        let gnomad = &mappings["gnomad"];
        assert_eq!(gnomad[0].src, "AF");
        assert_eq!(gnomad[0].dest, "gnomAD_AF");
        assert_eq!(gnomad[1].src, "AC");
        assert_eq!(gnomad[1].dest, "gnomAD_AC");
        Ok(())
    }

    #[test]
    fn annotate_vcf_adds_exact_key_info_fields() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tDP=5\nchr1\t11\t.\tA\tG\t.\tPASS\t.\n",
        )?;
        std::fs::write(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25;AC=3\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings =
            parse_annotation_mappings(&["db:AF=db_AF,AC=db_AC".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        annotate_vcf(
            &input,
            &output,
            IndexType::Csi,
            &databases,
            &mappings,
            &mut lookup,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;

        assert!(text.contains("##INFO=<ID=db_AF"));
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tDP=5;db_AC=3;db_AF=0.25"));
        assert!(text.contains("chr1\t11\t.\tA\tG\t.\tPASS\t."));
        assert!(dir.path().join("out.vcf.gz.csi").exists());
        Ok(())
    }

    #[test]
    fn annotate_vcf_writes_tbi_index_when_requested() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\t.\n",
        )?;
        std::fs::write(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        annotate_vcf(
            &input,
            &output,
            IndexType::Tbi,
            &databases,
            &mappings,
            &mut lookup,
        )?;

        assert!(dir.path().join("out.vcf.gz.tbi").exists());
        Ok(())
    }

    #[test]
    fn annotate_vcf_uses_tabix_indexed_database_when_available() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf.gz");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\t.\n",
        )?;
        write_bgzipped_vcf_with_tbi(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.5\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            &ExecutionContext {
                verbose: 0,
                threads: 1,
            },
        )?;
        annotate_vcf(
            &input,
            &output,
            IndexType::Csi,
            &databases,
            &mappings,
            &mut lookup,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tdb_AF=0.25"));
        Ok(())
    }

    #[test]
    fn skips_multi_alt_records_for_initial_exact_key_mode() -> Result<()> {
        let mappings = vec![FieldMapping {
            src: "AF".to_string(),
            dest: "db_AF".to_string(),
        }];
        let parsed = parse_database_record("chr1\t10\t.\tA\tC,G\t.\tPASS\tAF=0.1,0.2", &mappings)?;
        assert!(parsed.is_none());
        Ok(())
    }

    #[test]
    fn annotates_flag_info_as_one() -> Result<()> {
        let mappings = vec![FieldMapping {
            src: "COMMON".to_string(),
            dest: "db_COMMON".to_string(),
        }];
        let (_, values) =
            parse_database_record("chr1\t10\t.\tA\tC\t.\tPASS\tCOMMON", &mappings)?.unwrap();
        assert_eq!(values["db_COMMON"], "1");
        Ok(())
    }

    fn write_bgzipped_vcf_with_tbi(path: &Path, text: &str) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = bgzf::io::writer::Builder::default().build_from_writer(file);
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

        for line in text.lines() {
            if line.starts_with('#') {
                writeln!(writer, "{line}")?;
            } else {
                let record = index_record_from_line(line)?.context("missing record")?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{line}")?;
                let chunk_end = writer.virtual_position();
                indexer.add_record(
                    &record.chrom,
                    record.position,
                    record.position,
                    Chunk::new(chunk_start, chunk_end),
                )?;
            }
        }
        writer.try_finish()?;

        let index = indexer.build();
        let mut index_writer = tabix::io::Writer::new(File::create(tabix_index_path(path))?);
        index_writer.write_index(&index)?;
        Ok(())
    }
}
