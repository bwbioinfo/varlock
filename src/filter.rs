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

#[derive(Clone, Debug)]
struct SampleGtFilter {
    sample: String,
    gt: String,
}

#[derive(Clone, Debug)]
struct SampleMinDpFilter {
    sample: String,
    min_dp: u32,
}

#[derive(Clone, Debug)]
struct GroupGtFilter {
    group: String,
    gt: String,
}

#[derive(Clone, Debug)]
struct GroupMinDpFilter {
    group: String,
    min_dp: u32,
}

#[derive(Debug)]
struct FilterSpec {
    require_info: HashSet<String>,
    exclude_info: HashSet<String>,
    max_info: Vec<MaxInfoFilter>,
    sample_groups: HashMap<String, Vec<String>>,
    sample_has_alt: Vec<String>,
    sample_gt: Vec<SampleGtFilter>,
    sample_min_dp: Vec<SampleMinDpFilter>,
    group_any_has_alt: Vec<String>,
    group_any_gt: Vec<GroupGtFilter>,
    group_all_min_dp: Vec<GroupMinDpFilter>,
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
        let sample_groups = args
            .sample_groups
            .iter()
            .map(|value| {
                let (name, samples) = split_name_value(value, "--sample-group")?;
                let samples = samples
                    .split(',')
                    .filter(|sample| !sample.is_empty())
                    .map(str::to_string)
                    .collect::<Vec<_>>();
                if samples.is_empty() {
                    bail!("invalid --sample-group {value:?}; expected NAME=SAMPLE[,SAMPLE...]");
                }
                Ok((name.to_string(), samples))
            })
            .collect::<Result<HashMap<_, _>>>()?;
        let sample_gt = args
            .sample_gt
            .iter()
            .map(|value| {
                let (sample, gt) = split_name_value(value, "--sample-gt")?;
                Ok(SampleGtFilter {
                    sample: sample.to_string(),
                    gt: gt.to_string(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let sample_min_dp = args
            .sample_min_dp
            .iter()
            .map(|value| {
                let (sample, min_dp) = split_name_value(value, "--sample-min-dp")?;
                Ok(SampleMinDpFilter {
                    sample: sample.to_string(),
                    min_dp: parse_dp(min_dp, "--sample-min-dp")?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let group_any_gt = args
            .group_any_gt
            .iter()
            .map(|value| {
                let (group, gt) = split_name_value(value, "--group-any-gt")?;
                Ok(GroupGtFilter {
                    group: group.to_string(),
                    gt: gt.to_string(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let group_all_min_dp = args
            .group_all_min_dp
            .iter()
            .map(|value| {
                let (group, min_dp) = split_name_value(value, "--group-all-min-dp")?;
                Ok(GroupMinDpFilter {
                    group: group.to_string(),
                    min_dp: parse_dp(min_dp, "--group-all-min-dp")?,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            require_info: args.require_info.iter().cloned().collect(),
            exclude_info: args.exclude_info.iter().cloned().collect(),
            max_info,
            sample_groups,
            sample_has_alt: args.sample_has_alt.clone(),
            sample_gt,
            sample_min_dp,
            group_any_has_alt: args.group_any_has_alt.clone(),
            group_any_gt,
            group_all_min_dp,
        })
    }

    fn validate_samples(&self, sample_names: &[String]) -> Result<()> {
        let known = sample_names
            .iter()
            .map(String::as_str)
            .collect::<HashSet<_>>();
        for sample in &self.sample_has_alt {
            validate_known_sample(sample, &known)?;
        }
        for filter in &self.sample_gt {
            validate_known_sample(&filter.sample, &known)?;
        }
        for filter in &self.sample_min_dp {
            validate_known_sample(&filter.sample, &known)?;
        }
        for (group, samples) in &self.sample_groups {
            if samples.is_empty() {
                bail!("sample group {group:?} has no samples");
            }
            for sample in samples {
                validate_known_sample(sample, &known)
                    .with_context(|| format!("invalid sample group {group:?}"))?;
            }
        }
        for group in &self.group_any_has_alt {
            self.group_samples(group)?;
        }
        for filter in &self.group_any_gt {
            self.group_samples(&filter.group)?;
        }
        for filter in &self.group_all_min_dp {
            self.group_samples(&filter.group)?;
        }
        Ok(())
    }

    fn keep_record(&self, record: &VcfRecord<'_>) -> Result<bool> {
        let info = parse_info(record.info_text);
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
        for sample in &self.sample_has_alt {
            if !record.sample_has_alt(sample)? {
                return Ok(false);
            }
        }
        for filter in &self.sample_gt {
            if record.sample_field(&filter.sample, "GT") != Some(filter.gt.as_str()) {
                return Ok(false);
            }
        }
        for filter in &self.sample_min_dp {
            if !record.sample_min_dp(&filter.sample, filter.min_dp)? {
                return Ok(false);
            }
        }
        for group in &self.group_any_has_alt {
            let mut any = false;
            for sample in self.group_samples(group)? {
                if record.sample_has_alt(sample)? {
                    any = true;
                    break;
                }
            }
            if !any {
                return Ok(false);
            }
        }
        for filter in &self.group_any_gt {
            if !self
                .group_samples(&filter.group)?
                .iter()
                .any(|sample| record.sample_field(sample, "GT") == Some(filter.gt.as_str()))
            {
                return Ok(false);
            }
        }
        for filter in &self.group_all_min_dp {
            for sample in self.group_samples(&filter.group)? {
                if !record.sample_min_dp(sample, filter.min_dp)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn group_samples(&self, group: &str) -> Result<&[String]> {
        self.sample_groups
            .get(group)
            .map(Vec::as_slice)
            .with_context(|| format!("sample group {group:?} is not defined"))
    }
}

#[derive(Default)]
struct FilterMetrics {
    input_records: usize,
    output_records: usize,
}

struct VcfRecord<'a> {
    info_text: &'a str,
    fields: &'a [&'a str],
    sample_index: &'a HashMap<String, usize>,
    format_index: HashMap<&'a str, usize>,
}

impl<'a> VcfRecord<'a> {
    fn new(fields: &'a [&'a str], sample_index: &'a HashMap<String, usize>) -> Result<Self> {
        let format_index = if fields.len() >= 9 {
            fields[8]
                .split(':')
                .enumerate()
                .map(|(i, field)| (field, i))
                .collect()
        } else {
            HashMap::new()
        };
        Ok(Self {
            info_text: fields[7],
            fields,
            sample_index,
            format_index,
        })
    }

    fn sample_field(&self, sample: &str, field: &str) -> Option<&'a str> {
        let sample_offset = *self.sample_index.get(sample)?;
        let format_offset = *self.format_index.get(field)?;
        let sample_text = self.fields.get(9 + sample_offset)?;
        sample_text.split(':').nth(format_offset)
    }

    fn sample_has_alt(&self, sample: &str) -> Result<bool> {
        let Some(gt) = self.sample_field(sample, "GT") else {
            return Ok(false);
        };
        genotype_has_alt(gt)
    }

    fn sample_min_dp(&self, sample: &str, min_dp: u32) -> Result<bool> {
        let Some(dp) = self.sample_field(sample, "DP") else {
            return Ok(false);
        };
        let dp = dp
            .parse::<u32>()
            .with_context(|| format!("FORMAT/DP for sample {sample:?} is not an integer"))?;
        Ok(dp >= min_dp)
    }
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
    let mut sample_names = Vec::new();
    let mut sample_index = HashMap::new();
    let mut metrics = FilterMetrics::default();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            saw_column_header = true;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() > 9 {
                sample_names = fields[9..]
                    .iter()
                    .map(|sample| (*sample).to_string())
                    .collect();
                sample_index = sample_names
                    .iter()
                    .enumerate()
                    .map(|(i, sample)| (sample.clone(), i))
                    .collect();
            }
            spec.validate_samples(&sample_names)?;
            writeln!(writer, "{trimmed}")?;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.input_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let record = VcfRecord::new(&fields, &sample_index)?;
            if spec.keep_record(&record)? {
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

fn split_name_value<'a>(value: &'a str, flag: &str) -> Result<(&'a str, &'a str)> {
    let (name, parsed_value) = value
        .split_once('=')
        .with_context(|| format!("invalid {flag} {value:?}; expected NAME=VALUE"))?;
    if name.is_empty() || parsed_value.is_empty() {
        bail!("invalid {flag} {value:?}; expected NAME=VALUE");
    }
    Ok((name, parsed_value))
}

fn parse_dp(value: &str, flag: &str) -> Result<u32> {
    value
        .parse::<u32>()
        .with_context(|| format!("invalid {flag} depth threshold {value:?}"))
}

fn validate_known_sample(sample: &str, known: &HashSet<&str>) -> Result<()> {
    if known.contains(sample) {
        Ok(())
    } else {
        bail!("sample {sample:?} was not found in VCF header")
    }
}

fn genotype_has_alt(gt: &str) -> Result<bool> {
    if gt == "." || gt == "./." || gt == ".|." {
        return Ok(false);
    }
    for allele in gt.split(['/', '|']) {
        if allele == "." || allele.is_empty() {
            continue;
        }
        let allele = allele
            .parse::<u32>()
            .with_context(|| format!("FORMAT/GT allele {allele:?} is not an integer"))?;
        if allele > 0 {
            return Ok(true);
        }
    }
    Ok(false)
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

    fn empty_spec() -> FilterSpec {
        FilterSpec {
            require_info: HashSet::new(),
            exclude_info: HashSet::new(),
            max_info: Vec::new(),
            sample_groups: HashMap::new(),
            sample_has_alt: Vec::new(),
            sample_gt: Vec::new(),
            sample_min_dp: Vec::new(),
            group_any_has_alt: Vec::new(),
            group_any_gt: Vec::new(),
            group_all_min_dp: Vec::new(),
        }
    }

    fn record<'a>(
        fields: &'a [&'a str],
        sample_index: &'a HashMap<String, usize>,
    ) -> VcfRecord<'a> {
        VcfRecord::new(fields, sample_index).unwrap()
    }

    fn sample_index(samples: &[&str]) -> HashMap<String, usize> {
        samples
            .iter()
            .enumerate()
            .map(|(i, sample)| ((*sample).to_string(), i))
            .collect()
    }

    #[test]
    fn filter_spec_requires_and_excludes_info() -> Result<()> {
        let spec = FilterSpec {
            require_info: ["AF".to_string()].into_iter().collect(),
            exclude_info: ["COMMON".to_string()].into_iter().collect(),
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.1"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "DP=10"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.1;COMMON"],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_max_info_to_all_numeric_values() -> Result<()> {
        let spec = FilterSpec {
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.01,0.1"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.01,0.2"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "DP=10"],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_sample_gt_and_dp_predicates() -> Result<()> {
        let spec = FilterSpec {
            sample_has_alt: vec!["Tumor".to_string()],
            sample_gt: vec![SampleGtFilter {
                sample: "Normal".to_string(),
                gt: "0/0".to_string(),
            }],
            sample_min_dp: vec![SampleMinDpFilter {
                sample: "Tumor".to_string(),
                min_dp: 10,
            }],
            ..empty_spec()
        };
        let sample_index = sample_index(&["Tumor", "Normal"]);
        assert!(spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/1:12", "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/1:8", "0/0:20",
            ],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_group_predicates() -> Result<()> {
        let spec = FilterSpec {
            sample_groups: [
                (
                    "affected".to_string(),
                    vec!["TumorA".to_string(), "TumorB".to_string()],
                ),
                ("controls".to_string(), vec!["Normal".to_string()]),
            ]
            .into_iter()
            .collect(),
            group_any_has_alt: vec!["affected".to_string()],
            group_any_gt: vec![GroupGtFilter {
                group: "affected".to_string(),
                gt: "0/1".to_string(),
            }],
            group_all_min_dp: vec![GroupMinDpFilter {
                group: "controls".to_string(),
                min_dp: 15,
            }],
            ..empty_spec()
        };
        let sample_index = sample_index(&["TumorA", "TumorB", "Normal"]);
        assert!(spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/1:14",
                "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/0:14",
                "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/1:14",
                "0/0:10",
            ],
            &sample_index
        ))?);
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
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
            ..empty_spec()
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

    #[test]
    fn filter_vcf_applies_sample_aware_predicates() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let output = dir.path().join("out.vcf.gz");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal\n\
             chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05\tGT:DP\t0/1:12\t0/0:20\n\
             chr1\t11\t.\tA\tG\t.\tPASS\tAF=0.05\tGT:DP\t0/0:12\t0/0:20\n\
             chr1\t12\t.\tA\tT\t.\tPASS\tAF=0.05\tGT:DP\t0/1:8\t0/0:20\n",
        )?;
        let spec = FilterSpec {
            sample_has_alt: vec!["Tumor".to_string()],
            sample_min_dp: vec![SampleMinDpFilter {
                sample: "Tumor".to_string(),
                min_dp: 10,
            }],
            sample_gt: vec![SampleGtFilter {
                sample: "Normal".to_string(),
                gt: "0/0".to_string(),
            }],
            ..empty_spec()
        };

        let metrics = filter_vcf(&input, &output, IndexType::Csi, &spec)?;
        assert_eq!(metrics.input_records, 3);
        assert_eq!(metrics.output_records, 1);

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(!text.contains("chr1\t12\t.\tA\tT"));
        Ok(())
    }
}
