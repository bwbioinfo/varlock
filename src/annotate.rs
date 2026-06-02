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
use noodles_csi::{
    self as csi,
    binning_index::{
        self,
        index::reference_sequence::bin::Chunk,
        index::reference_sequence::index::BinnedIndex,
        index::{
            Header as TabixHeader,
            header::{Format as TabixFormat, ReferenceSequenceNames},
        },
    },
};
use noodles_tabix as tabix;

use crate::call_targets::reference::{FastaIndex, open_fasta_index};
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

struct VariantNormalizer {
    fasta: FastaIndex,
}

impl VariantNormalizer {
    fn open(reference: &Path) -> Result<Self> {
        Ok(Self {
            fasta: open_fasta_index(reference)?,
        })
    }

    fn normalize_key(&mut self, key: VariantKey) -> Result<VariantKey> {
        if is_symbolic_or_special_allele(&key.ref_allele)
            || is_symbolic_or_special_allele(&key.alt_allele)
        {
            return Ok(key);
        }

        let mut pos = key
            .pos
            .parse::<u32>()
            .with_context(|| format!("invalid VCF position {:?}", key.pos))?;
        let mut ref_allele = key.ref_allele.to_ascii_uppercase();
        let mut alt_allele = key.alt_allele.to_ascii_uppercase();

        trim_common_suffix(&mut ref_allele, &mut alt_allele);
        trim_common_prefix(&mut pos, &mut ref_allele, &mut alt_allele);

        if ref_allele.len() != alt_allele.len() {
            while pos > 1 {
                let prev_base = self.fasta.fetch_base(&key.chrom, pos - 1)? as char;
                let Some(ref_last) = ref_allele.chars().last() else {
                    break;
                };
                let Some(alt_last) = alt_allele.chars().last() else {
                    break;
                };
                if ref_last != prev_base || alt_last != prev_base {
                    break;
                }
                ref_allele.pop();
                alt_allele.pop();
                ref_allele.insert(0, prev_base);
                alt_allele.insert(0, prev_base);
                pos -= 1;
            }
        }

        trim_common_suffix(&mut ref_allele, &mut alt_allele);
        trim_common_prefix(&mut pos, &mut ref_allele, &mut alt_allele);

        Ok(VariantKey {
            chrom: key.chrom,
            pos: pos.to_string(),
            ref_allele,
            alt_allele,
        })
    }
}

fn is_symbolic_or_special_allele(allele: &str) -> bool {
    allele == "*"
        || allele.starts_with('<')
        || allele.contains('>')
        || allele.contains('[')
        || allele.contains(']')
}

fn trim_common_suffix(ref_allele: &mut String, alt_allele: &mut String) {
    while ref_allele.len() > 1 && alt_allele.len() > 1 {
        let Some(ref_last) = ref_allele.chars().last() else {
            break;
        };
        let Some(alt_last) = alt_allele.chars().last() else {
            break;
        };
        if ref_last != alt_last {
            break;
        }
        ref_allele.pop();
        alt_allele.pop();
    }
}

fn trim_common_prefix(pos: &mut u32, ref_allele: &mut String, alt_allele: &mut String) {
    while ref_allele.len() > 1 && alt_allele.len() > 1 {
        let Some(ref_first) = ref_allele.chars().next() else {
            break;
        };
        let Some(alt_first) = alt_allele.chars().next() else {
            break;
        };
        if ref_first != alt_first {
            break;
        }
        ref_allele.remove(0);
        alt_allele.remove(0);
        *pos = pos.saturating_add(1);
    }
}

pub(crate) fn run(args: AnnotateArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    let databases = parse_database_specs(&args.databases)?;
    let mappings = parse_annotation_mappings(&args.annotations, &databases)?;
    let mut normalizer = args
        .reference
        .as_deref()
        .map(VariantNormalizer::open)
        .transpose()?;
    let mut lookup = AnnotationLookup::open(&databases, &mappings, normalizer.as_mut(), ctx)?;

    annotate_vcf(
        &args.input,
        &args.output,
        args.index_type,
        &databases,
        &mappings,
        &mut lookup,
        normalizer.as_mut(),
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
    Csi {
        name: String,
        mappings: Vec<FieldMapping>,
        reader: noodles_csi::io::IndexedReader<bgzf::io::Reader<File>, csi::Index>,
    },
}

struct AnnotationLookup {
    databases: Vec<DatabaseLookup>,
}

impl AnnotationLookup {
    fn open(
        databases: &[DatabaseSpec],
        mappings: &HashMap<String, Vec<FieldMapping>>,
        mut normalizer: Option<&mut VariantNormalizer>,
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

            if normalizer.is_none() && tabix_index_path(&db.path).exists() {
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
            } else if normalizer.is_none() && csi_index_path(&db.path).exists() {
                let index_path = csi_index_path(&db.path);
                let index = csi::fs::read(&index_path).with_context(|| {
                    format!("failed to read CSI index {}", index_path.display())
                })?;
                let file = File::open(&db.path).with_context(|| {
                    format!(
                        "failed to open CSI-indexed annotation database {}",
                        db.path.display()
                    )
                })?;
                let reader = csi::io::IndexedReader::new(file, index);
                log_verbose(
                    ctx,
                    format!(
                        "annotate stage=open_database name={} mode=csi path={}",
                        db.name,
                        db.path.display()
                    ),
                );
                opened.push(DatabaseLookup::Csi {
                    name: db.name.clone(),
                    mappings: db_mappings.clone(),
                    reader,
                });
            } else {
                let annotations =
                    load_one_annotation_database(db, db_mappings, normalizer.as_deref_mut())?;
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
                    out.extend(query_indexed_database(name, mappings, reader, key)?);
                }
                DatabaseLookup::Csi {
                    name,
                    mappings,
                    reader,
                } => {
                    out.extend(query_indexed_database(name, mappings, reader, key)?);
                }
            }
        }
        Ok(out)
    }

    fn lookup_record(&mut self, keys: &[VariantKey]) -> Result<AnnotationValues> {
        let mut per_alt = Vec::with_capacity(keys.len());
        for key in keys {
            per_alt.push(self.lookup(key)?);
        }

        let mut field_names = BTreeMap::new();
        for values in &per_alt {
            for field in values.keys() {
                field_names.insert(field.clone(), ());
            }
        }

        let mut out = AnnotationValues::new();
        for field in field_names.keys() {
            let values = per_alt
                .iter()
                .map(|values| {
                    values
                        .get(field)
                        .cloned()
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect::<Vec<_>>();
            if values.iter().any(|value| value != ".") {
                out.insert(field.clone(), values.join(","));
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
                DatabaseLookup::Csi { .. } => 0,
            })
            .sum()
    }
}

fn load_one_annotation_database(
    db: &DatabaseSpec,
    db_mappings: &[FieldMapping],
    mut normalizer: Option<&mut VariantNormalizer>,
) -> Result<AnnotationMap> {
    let mut out: AnnotationMap = HashMap::new();
    let mut reader = open_text_reader(&db.path)
        .with_context(|| format!("failed to open annotation database {}", db.path.display()))?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            let records = parse_database_records(trimmed, db_mappings, normalizer.as_deref_mut())?;
            for (key, values) in records {
                out.entry(key).or_default().extend(values);
            }
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
    mut normalizer: Option<&mut VariantNormalizer>,
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
            let annotated = annotate_record_line(trimmed, lookup, normalizer.as_deref_mut())?;
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
                    "##INFO=<ID={},Number=A,Type=String,Description=\"Annotation from {}:{}\">",
                    mapping.dest, db.name, mapping.src
                )?;
            }
        }
    }
    Ok(())
}

fn annotate_record_line(
    line: &str,
    lookup: &mut AnnotationLookup,
    normalizer: Option<&mut VariantNormalizer>,
) -> Result<String> {
    let mut fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid VCF record with fewer than 8 fields: {line}");
    }
    let keys = variant_keys_from_fields(&fields, normalizer)?;
    if keys.is_empty() {
        return Ok(line.to_string());
    }
    let values = lookup.lookup_record(&keys)?;
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

fn query_indexed_database<I>(
    db_name: &str,
    mappings: &[FieldMapping],
    reader: &mut noodles_csi::io::IndexedReader<bgzf::io::Reader<File>, I>,
    key: &VariantKey,
) -> Result<AnnotationValues>
where
    I: csi::BinningIndex,
{
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
        for (record_key, values) in parse_database_records(line, mappings, None)? {
            if record_key == *key {
                out.extend(values);
            }
        }
    }
    Ok(out)
}

fn parse_database_records(
    line: &str,
    mappings: &[FieldMapping],
    normalizer: Option<&mut VariantNormalizer>,
) -> Result<Vec<(VariantKey, AnnotationValues)>> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid database VCF record with fewer than 8 fields: {line}");
    }
    let keys = variant_keys_from_fields(&fields, normalizer)?;
    if keys.is_empty() {
        return Ok(Vec::new());
    }
    let alt_count = keys.len();
    let info = parse_info(fields[7]);
    let mut records = Vec::new();
    for (alt_index, key) in keys.into_iter().enumerate() {
        let mut values = AnnotationValues::new();
        for mapping in mappings {
            if let Some(value) = info.get(mapping.src.as_str()) {
                values.insert(
                    mapping.dest.clone(),
                    annotation_value_for_alt(value, alt_count, alt_index),
                );
            }
        }
        if !values.is_empty() {
            records.push((key, values));
        }
    }
    Ok(records)
}

fn variant_keys_from_fields(
    fields: &[&str],
    mut normalizer: Option<&mut VariantNormalizer>,
) -> Result<Vec<VariantKey>> {
    let Some(chrom) = fields.first() else {
        return Ok(Vec::new());
    };
    let Some(pos) = fields.get(1) else {
        return Ok(Vec::new());
    };
    let Some(ref_allele) = fields.get(3) else {
        return Ok(Vec::new());
    };
    let Some(alt_field) = fields.get(4) else {
        return Ok(Vec::new());
    };
    if *alt_field == "." || alt_field.is_empty() {
        return Ok(Vec::new());
    }

    let mut keys = Vec::new();
    for alt in alt_field
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
    {
        let key = VariantKey {
            chrom: (*chrom).to_string(),
            pos: (*pos).to_string(),
            ref_allele: (*ref_allele).to_string(),
            alt_allele: alt.to_string(),
        };
        let key = if let Some(normalizer) = normalizer.as_deref_mut() {
            normalizer.normalize_key(key)?
        } else {
            key
        };
        keys.push(key);
    }
    Ok(keys)
}

fn annotation_value_for_alt(value: &str, alt_count: usize, alt_index: usize) -> String {
    let parts = value.split(',').collect::<Vec<_>>();
    if parts.len() == alt_count {
        parts
            .get(alt_index)
            .copied()
            .unwrap_or_default()
            .to_string()
    } else {
        value.to_string()
    }
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

fn csi_index_path(path: &Path) -> PathBuf {
    path.with_file_name(format!(
        "{}.csi",
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
            None,
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
            None,
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
            None,
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
            None,
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
            None,
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
            None,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tdb_AF=0.25"));
        Ok(())
    }

    #[test]
    fn annotate_vcf_uses_csi_indexed_database_when_available() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf.gz");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\t.\n",
        )?;
        write_bgzipped_vcf_with_csi(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.25\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.5\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            None,
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
            None,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tdb_AF=0.25"));
        Ok(())
    }

    #[test]
    fn parses_multi_alt_database_records_by_alt() -> Result<()> {
        let mappings = vec![FieldMapping {
            src: "AF".to_string(),
            dest: "db_AF".to_string(),
        }];
        let parsed =
            parse_database_records("chr1\t10\t.\tA\tC,G\t.\tPASS\tAF=0.1,0.2", &mappings, None)?;
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].0.alt_allele, "C");
        assert_eq!(parsed[0].1["db_AF"], "0.1");
        assert_eq!(parsed[1].0.alt_allele, "G");
        assert_eq!(parsed[1].1["db_AF"], "0.2");
        Ok(())
    }

    #[test]
    fn annotate_vcf_adds_multi_alt_info_in_alt_order() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC,G\t.\tPASS\t.\n",
        )?;
        std::fs::write(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC,G\t.\tPASS\tAF=0.1,0.2\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            None,
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
            None,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("##INFO=<ID=db_AF,Number=A"));
        assert!(text.contains("chr1\t10\t.\tA\tC,G\t.\tPASS\tdb_AF=0.1,0.2"));
        Ok(())
    }

    #[test]
    fn annotate_vcf_marks_missing_multi_alt_annotations_with_dot() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf");
        let output = dir.path().join("out.vcf.gz");

        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC,G\t.\tPASS\t.\n",
        )?;
        std::fs::write(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.1\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            None,
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
            None,
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC,G\t.\tPASS\tdb_AF=0.1,."));
        Ok(())
    }

    #[test]
    fn annotate_vcf_matches_left_shifted_indel_with_reference() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let db = dir.path().join("db.vcf");
        let output = dir.path().join("out.vcf.gz");
        let reference = dir.path().join("ref.fa");

        write_reference_with_fai(&reference, "chr1", "AAAAAA")?;
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\t.\tAA\tA\t.\tPASS\t.\n",
        )?;
        std::fs::write(
            &db,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tAA\tA\t.\tPASS\tAF=0.25\n",
        )?;

        let databases = parse_database_specs(&[format!("db={}", db.display())])?;
        let mappings = parse_annotation_mappings(&["db:AF=db_AF".to_string()], &databases)?;
        let mut normalizer = VariantNormalizer::open(&reference)?;
        let mut lookup = AnnotationLookup::open(
            &databases,
            &mappings,
            Some(&mut normalizer),
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
            Some(&mut normalizer),
        )?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t3\t.\tAA\tA\t.\tPASS\tdb_AF=0.25"));
        Ok(())
    }

    #[test]
    fn normalizer_leaves_symbolic_alleles_unchanged() -> Result<()> {
        let dir = tempdir()?;
        let reference = dir.path().join("ref.fa");
        write_reference_with_fai(&reference, "chr1", "AAAAAA")?;
        let mut normalizer = VariantNormalizer::open(&reference)?;
        let key = VariantKey {
            chrom: "chr1".to_string(),
            pos: "3".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "<DEL>".to_string(),
        };
        let normalized = normalizer.normalize_key(key.clone())?;
        assert_eq!(normalized, key);
        Ok(())
    }

    #[test]
    fn annotates_flag_info_as_one() -> Result<()> {
        let mappings = vec![FieldMapping {
            src: "COMMON".to_string(),
            dest: "db_COMMON".to_string(),
        }];
        let parsed = parse_database_records("chr1\t10\t.\tA\tC\t.\tPASS\tCOMMON", &mappings, None)?;
        let (_, values) = parsed.first().context("missing parsed record")?;
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

    fn write_bgzipped_vcf_with_csi(path: &Path, text: &str) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = bgzf::io::writer::Builder::default().build_from_writer(file);
        let mut indexer = binning_index::Indexer::<BinnedIndex>::default();
        let mut reference_ids = HashMap::new();
        let mut reference_names = Vec::new();

        for line in text.lines() {
            if line.starts_with('#') {
                writeln!(writer, "{line}")?;
            } else {
                let record = index_record_from_line(line)?.context("missing record")?;
                let reference_sequence_id =
                    reference_id_for(&mut reference_ids, &mut reference_names, &record.chrom);
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{line}")?;
                let chunk_end = writer.virtual_position();
                indexer.add_record(
                    Some((
                        reference_sequence_id,
                        record.position,
                        record.position,
                        true,
                    )),
                    Chunk::new(chunk_start, chunk_end),
                )?;
            }
        }
        writer.try_finish()?;

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
        let mut index_writer = csi::io::Writer::new(File::create(csi_index_path(path))?);
        index_writer.write_index(&index)?;
        Ok(())
    }

    fn write_reference_with_fai(path: &Path, name: &str, sequence: &str) -> Result<()> {
        let fasta = format!(">{name}\n{sequence}\n");
        std::fs::write(path, fasta)?;
        let offset = name.len() + 2;
        let fai = format!(
            "{name}\t{}\t{offset}\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        );
        std::fs::write(path.with_extension("fa.fai"), fai)?;
        Ok(())
    }
}
