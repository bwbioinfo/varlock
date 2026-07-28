use std::{
    collections::BTreeSet,
    fs::{self, File},
    io::ErrorKind,
    path::{Path, PathBuf},
};

#[cfg(unix)]
use std::os::unix::fs::MetadataExt;

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf, io::Write as _, record::data::field::Tag, record_buf::data::field::Value,
    },
    header::record::value::{
        Map,
        map::{ReadGroup, read_group::tag::SAMPLE},
    },
};

use crate::{AddSmToBamArgs, ExecutionContext, log_verbose};

pub(crate) fn run(args: AddSmToBamArgs, ctx: &ExecutionContext) -> Result<()> {
    if paths_refer_to_same_file(&args.input, &args.output)? {
        bail!("--output must differ from --input when rewriting a BAM header");
    }

    let sample = args
        .sample
        .unwrap_or_else(|| default_sample_name(&args.input));
    validate_sample_name(&sample)?;

    let input = File::open(&args.input)
        .with_context(|| format!("failed to open input BAM {}", args.input.display()))?;
    let mut reader = bam::io::Reader::new(input);
    let mut header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header {}", args.input.display()))?;
    let read_group_update = prepare_read_groups(&mut header, &args.input, &sample)?;

    let output = File::create(&args.output)
        .with_context(|| format!("failed to create output BAM {}", args.output.display()))?;
    let mut writer = bam::io::Writer::new(output);
    writer
        .write_header(&header)
        .with_context(|| format!("failed to write BAM header {}", args.output.display()))?;

    let mut record_count = 0u64;
    let mut tagged_record_count = 0u64;
    for result in reader.records() {
        let record = result
            .with_context(|| format!("failed to read record from {}", args.input.display()))?;

        if let Some(read_group_id) = read_group_update.default_read_group_for_untagged.as_deref() {
            let mut record =
                RecordBuf::try_from_alignment_record(&header, &record).with_context(|| {
                    format!(
                        "failed to buffer record from {} while adding its RG tag",
                        args.input.display()
                    )
                })?;

            if record_read_group_id(&record, &args.input)?.is_none() {
                record
                    .data_mut()
                    .insert(Tag::READ_GROUP, Value::from(read_group_id));
                tagged_record_count += 1;
            }

            writer
                .write_alignment_record(&header, &record)
                .with_context(|| format!("failed to write record to {}", args.output.display()))?;
        } else {
            writer
                .write_alignment_record(&header, &record)
                .with_context(|| format!("failed to write record to {}", args.output.display()))?;
        }

        record_count += 1;
    }
    writer
        .try_finish()
        .with_context(|| format!("failed to finish output BAM {}", args.output.display()))?;
    drop(writer);

    let index = bam::fs::index(&args.output)
        .with_context(|| format!("failed to index output BAM {}", args.output.display()))?;
    let index_path = bai_path(&args.output);
    bam::bai::fs::write(&index_path, &index)
        .with_context(|| format!("failed to write BAM index {}", index_path.display()))?;

    log_verbose(
        ctx,
        format!(
            "add-sm-to-bam input={} output={} sample={} sm_tags_added={} read_groups_added={} records_tagged={} records={}",
            args.input.display(),
            args.output.display(),
            sample,
            read_group_update.sample_tags_added,
            read_group_update.read_groups_added,
            tagged_record_count,
            record_count,
        ),
    );

    Ok(())
}

fn default_sample_name(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .filter(|stem| !stem.is_empty())
        .unwrap_or("sample")
        .to_owned()
}

fn default_read_group_id(path: &Path) -> String {
    default_sample_name(path)
}

fn validate_sample_name(sample: &str) -> Result<()> {
    if sample.is_empty() {
        bail!("--sample must not be empty");
    }
    if sample
        .bytes()
        .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r'))
    {
        bail!("--sample must not contain tabs or line breaks");
    }
    Ok(())
}

fn add_missing_sample_tags(header: &mut sam::Header, sample: &str) -> usize {
    let mut updated = 0usize;
    for read_group in header.read_groups_mut().values_mut() {
        if read_group.other_fields().contains_key(&SAMPLE) {
            continue;
        }
        read_group
            .other_fields_mut()
            .insert(SAMPLE, sample.to_owned().into());
        updated += 1;
    }

    updated
}

#[derive(Debug, Default)]
struct ReadGroupUpdate {
    sample_tags_added: usize,
    read_groups_added: usize,
    default_read_group_for_untagged: Option<String>,
}

#[derive(Debug, Default)]
struct RecordReadGroups {
    ids: BTreeSet<String>,
    has_untagged_records: bool,
}

fn prepare_read_groups(
    header: &mut sam::Header,
    input: &Path,
    sample: &str,
) -> Result<ReadGroupUpdate> {
    if !header.read_groups().is_empty() {
        return Ok(ReadGroupUpdate {
            sample_tags_added: add_missing_sample_tags(header, sample),
            ..ReadGroupUpdate::default()
        });
    }

    let mut record_read_groups = collect_record_read_groups(input)?;
    let default_read_group_id = default_read_group_id(input);
    let needs_default_read_group =
        record_read_groups.has_untagged_records || record_read_groups.ids.is_empty();

    if needs_default_read_group {
        record_read_groups.ids.insert(default_read_group_id.clone());
    }

    for read_group_id in &record_read_groups.ids {
        add_read_group(header, read_group_id, sample)?;
    }

    Ok(ReadGroupUpdate {
        sample_tags_added: record_read_groups.ids.len(),
        read_groups_added: record_read_groups.ids.len(),
        default_read_group_for_untagged: needs_default_read_group.then_some(default_read_group_id),
    })
}

fn collect_record_read_groups(path: &Path) -> Result<RecordReadGroups> {
    let input =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = bam::io::Reader::new(input);
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header {}", path.display()))?;
    let mut read_groups = RecordReadGroups::default();

    for result in reader.records() {
        let record =
            result.with_context(|| format!("failed to read record from {}", path.display()))?;
        let record = RecordBuf::try_from_alignment_record(&header, &record)
            .with_context(|| format!("failed to buffer record from {}", path.display()))?;

        match record_read_group_id(&record, path)? {
            Some(read_group_id) => {
                read_groups.ids.insert(read_group_id);
            }
            None => read_groups.has_untagged_records = true,
        }
    }

    Ok(read_groups)
}

fn add_read_group(header: &mut sam::Header, read_group_id: &str, sample: &str) -> Result<()> {
    validate_read_group_id(read_group_id)?;
    let read_group = Map::<ReadGroup>::builder()
        .insert(SAMPLE, sample)
        .build()
        .context("failed to build synthesized @RG record")?;
    header
        .read_groups_mut()
        .insert(read_group_id.into(), read_group);
    Ok(())
}

fn record_read_group_id(record: &RecordBuf, source: &Path) -> Result<Option<String>> {
    let value = match record.data().get(&Tag::READ_GROUP) {
        Some(Value::String(value)) | Some(Value::Hex(value)) => value,
        Some(_) => bail!("RG tag has unexpected type in {}", source.display()),
        None => return Ok(None),
    };
    let read_group_id = std::str::from_utf8(value.as_ref())
        .with_context(|| format!("invalid RG tag in {}", source.display()))?;
    validate_read_group_id(read_group_id)?;
    Ok(Some(read_group_id.to_owned()))
}

fn validate_read_group_id(read_group_id: &str) -> Result<()> {
    if read_group_id.is_empty() {
        bail!("read-group ID must not be empty");
    }
    if read_group_id
        .bytes()
        .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r'))
    {
        bail!("read-group ID must not contain tabs or line breaks");
    }
    Ok(())
}

fn bai_path(path: &Path) -> PathBuf {
    let mut index_path = path.as_os_str().to_os_string();
    index_path.push(".bai");
    PathBuf::from(index_path)
}

fn paths_refer_to_same_file(input: &Path, output: &Path) -> Result<bool> {
    let input_canonical = fs::canonicalize(input)
        .with_context(|| format!("failed to resolve input BAM {}", input.display()))?;
    let output_canonical = match fs::canonicalize(output) {
        Ok(path) => path,
        Err(error) if error.kind() == ErrorKind::NotFound => return Ok(false),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("failed to resolve output BAM {}", output.display()));
        }
    };

    if input_canonical == output_canonical {
        return Ok(true);
    }

    #[cfg(unix)]
    {
        let input_metadata = fs::metadata(input)
            .with_context(|| format!("failed to inspect input BAM {}", input.display()))?;
        let output_metadata = fs::metadata(output)
            .with_context(|| format!("failed to inspect output BAM {}", output.display()))?;
        Ok(input_metadata.dev() == output_metadata.dev()
            && input_metadata.ino() == output_metadata.ino())
    }

    #[cfg(not(unix))]
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::{
        add_missing_sample_tags, bai_path, default_read_group_id, default_sample_name,
        paths_refer_to_same_file, validate_read_group_id, validate_sample_name,
    };
    use anyhow::Result;
    use noodles_sam::{
        self as sam,
        header::record::value::{Map, map::ReadGroup},
    };
    use std::{
        fs::File,
        path::{Path, PathBuf},
    };
    use tempfile::tempdir;

    #[test]
    fn default_sample_name_uses_input_filename_stem() {
        assert_eq!(default_sample_name(Path::new("/data/run-01.bam")), "run-01");
        assert_eq!(
            default_read_group_id(Path::new("/data/run-01.bam")),
            "run-01"
        );
        assert_eq!(default_sample_name(Path::new("/")), "sample");
    }

    #[test]
    fn add_missing_sample_tags_preserves_existing_values() -> Result<()> {
        use noodles_sam::header::record::value::map::read_group::tag::SAMPLE;

        let mut header = sam::Header::builder()
            .add_read_group("missing", Map::<ReadGroup>::default())
            .add_read_group(
                "existing",
                Map::<ReadGroup>::builder()
                    .insert(SAMPLE, "already-set")
                    .build()?,
            )
            .build();

        assert_eq!(add_missing_sample_tags(&mut header, "derived"), 1);
        assert_eq!(
            header.read_groups()[&b"missing"[..]]
                .other_fields()
                .get(&SAMPLE)
                .map(|value| value.as_ref()),
            Some(&b"derived"[..])
        );
        assert_eq!(
            header.read_groups()[&b"existing"[..]]
                .other_fields()
                .get(&SAMPLE)
                .map(|value| value.as_ref()),
            Some(&b"already-set"[..])
        );

        Ok(())
    }

    #[test]
    fn add_missing_sample_tags_leaves_empty_headers_unchanged() {
        let mut header = sam::Header::default();
        assert_eq!(add_missing_sample_tags(&mut header, "sample"), 0);
    }

    #[test]
    fn sample_validation_rejects_empty_and_header_delimiters() {
        assert!(validate_sample_name("").is_err());
        assert!(validate_sample_name("bad\tvalue").is_err());
        assert!(validate_sample_name("valid sample").is_ok());
    }

    #[test]
    fn read_group_validation_rejects_empty_and_header_delimiters() {
        assert!(validate_read_group_id("").is_err());
        assert!(validate_read_group_id("bad\nvalue").is_err());
        assert!(validate_read_group_id("lane-1").is_ok());
    }

    #[test]
    fn bai_path_appends_to_bam_extension() {
        assert_eq!(
            bai_path(Path::new("calls.bam")),
            PathBuf::from("calls.bam.bai")
        );
    }

    #[test]
    fn paths_refer_to_same_file_detects_path_aliases() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.bam");
        File::create(&input)?;
        let alias = dir.path().join(".").join("input.bam");

        assert!(paths_refer_to_same_file(&input, &alias)?);
        assert!(!paths_refer_to_same_file(
            &input,
            &dir.path().join("output.bam")
        )?);
        Ok(())
    }
}
