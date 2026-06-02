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

use crate::{AnnotateArgs, ExecutionContext, log_verbose};

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
    let annotations = load_annotation_databases(&databases, &mappings)?;

    annotate_vcf(
        &args.input,
        &args.output,
        &databases,
        &mappings,
        &annotations,
    )?;

    log_verbose(
        ctx,
        format!(
            "annotate stage=done databases={} annotated_sites={} output={} elapsed={:.2?}",
            databases.len(),
            annotations.len(),
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

fn load_annotation_databases(
    databases: &[DatabaseSpec],
    mappings: &HashMap<String, Vec<FieldMapping>>,
) -> Result<AnnotationMap> {
    let mut out: AnnotationMap = HashMap::new();
    for db in databases {
        let db_mappings = mappings.get(&db.name).with_context(|| {
            format!(
                "no --annotation mapping provided for database {:?}",
                db.name
            )
        })?;
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
    }
    Ok(out)
}

fn annotate_vcf(
    input: &Path,
    output: &Path,
    databases: &[DatabaseSpec],
    mappings: &HashMap<String, Vec<FieldMapping>>,
    annotations: &AnnotationMap,
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
            let annotated = annotate_record_line(trimmed, annotations)?;
            writeln!(writer, "{annotated}")?;
        }
        line.clear();
    }

    if !wrote_annotation_headers {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
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

fn annotate_record_line(line: &str, annotations: &AnnotationMap) -> Result<String> {
    let mut fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 8 {
        bail!("invalid VCF record with fewer than 8 fields: {line}");
    }
    let Some(key) = variant_key_from_fields(&fields) else {
        return Ok(line.to_string());
    };
    let Some(values) = annotations.get(&key) else {
        return Ok(line.to_string());
    };
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
        info.push_str(key);
        info.push('=');
        info.push_str(value);
    }
    fields[7] = &info;
    Ok(fields.join("\t"))
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

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use std::io::Read;
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
        let annotations = load_annotation_databases(&databases, &mappings)?;
        annotate_vcf(&input, &output, &databases, &mappings, &annotations)?;

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;

        assert!(text.contains("##INFO=<ID=db_AF"));
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tDP=5;db_AC=3;db_AF=0.25"));
        assert!(text.contains("chr1\t11\t.\tA\tG\t.\tPASS\t."));
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
}
