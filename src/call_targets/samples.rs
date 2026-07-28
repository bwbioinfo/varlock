use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use noodles_bam::{Record, io::Reader};
use noodles_sam::{
    self as sam,
    alignment::record::data::field::{Tag, Value},
};

use crate::SplitBy;

pub(crate) const DEFAULT_SAMPLE_NAME: &str = "SAMPLE";

#[derive(Clone, Debug)]
pub(crate) struct SampleResolution {
    pub(crate) sample_names: Vec<String>,
    pub(crate) sample_index: HashMap<String, usize>,
    pub(crate) input_resolvers: Vec<InputSampleResolver>,
}

#[derive(Clone, Debug)]
pub(crate) struct InputSampleResolver {
    source: PathBuf,
    fallback_sample_index: Option<usize>,
    rg_to_sample_index: HashMap<String, usize>,
    require_read_group: bool,
}

impl InputSampleResolver {
    pub(crate) fn resolve_record(&self, record: &Record) -> Result<Option<usize>> {
        if self.rg_to_sample_index.is_empty() {
            return Ok(self.fallback_sample_index);
        }

        let read_group = match record.data().get(&Tag::READ_GROUP) {
            None if self.require_read_group => {
                bail!(
                    "--split-by rg requires an RG tag on retained reads in {}",
                    self.source.display()
                );
            }
            None => return Ok(self.fallback_sample_index),
            Some(Ok(Value::String(value))) | Some(Ok(Value::Hex(value))) => {
                std::str::from_utf8(value.as_ref()).context("invalid RG tag")?
            }
            Some(Ok(_)) => bail!("RG tag has unexpected type in {}", self.source.display()),
            Some(Err(e)) => {
                return Err(e).with_context(|| {
                    format!("failed to read RG tag in {}", self.source.display())
                });
            }
        };

        match self.rg_to_sample_index.get(read_group).copied() {
            Some(sample_index) => Ok(Some(sample_index)),
            None if self.require_read_group => bail!(
                "--split-by rg found RG {:?} in {} that is not declared in its BAM header",
                read_group,
                self.source.display()
            ),
            None => Ok(self.fallback_sample_index),
        }
    }

    pub(crate) fn describe(&self, sample_names: &[String]) -> String {
        let mut mappings = self
            .rg_to_sample_index
            .iter()
            .map(|(rg, sample_index)| format!("{rg}->{}", sample_names[*sample_index]))
            .collect::<Vec<_>>();
        mappings.sort();
        if let Some(sample_index) = self.fallback_sample_index {
            mappings.push(format!("fallback->{}", sample_names[sample_index]));
        }
        if mappings.is_empty() {
            mappings.push("no resolved samples".to_string());
        }
        format!("{} -> {}", self.source.display(), mappings.join(","))
    }
}

pub(crate) fn collect_sample_resolution(
    paths: &[PathBuf],
    split_by: SplitBy,
    rg_map: Option<&[(String, String)]>,
) -> Result<SampleResolution> {
    if rg_map.is_some() && split_by != SplitBy::Sm {
        bail!("--rg-map is only supported with --split-by sm");
    }

    match split_by {
        SplitBy::Sm => collect_sm_resolution(paths, rg_map),
        SplitBy::File => collect_file_resolution(paths),
        SplitBy::Rg => collect_rg_resolution(paths),
    }
}

fn collect_sm_resolution(
    paths: &[PathBuf],
    rg_map: Option<&[(String, String)]>,
) -> Result<SampleResolution> {
    if let Some(rows) = rg_map {
        let (sample_names, rg_to_sm) = build_rg_to_sm_from_map(rows)?;
        let sample_index = sample_index(&sample_names);
        let rg_to_sample_index = rg_to_sm
            .into_iter()
            .map(|(rg, sm)| (rg, sample_index[&sm]))
            .collect::<HashMap<_, _>>();
        let input_resolvers = paths
            .iter()
            .cloned()
            .map(|source| InputSampleResolver {
                source,
                fallback_sample_index: None,
                rg_to_sample_index: rg_to_sample_index.clone(),
                require_read_group: false,
            })
            .collect();
        return Ok(SampleResolution {
            sample_names,
            sample_index,
            input_resolvers,
        });
    }

    let input_mappings = paths
        .iter()
        .map(|path| {
            let header = read_header(path)?;
            let rg_to_sm = rg_to_sm_from_header(&header)?;
            let read_group_ids = read_group_ids_from_header(&header)?;
            let has_unmapped_read_groups = read_group_ids.is_empty()
                || read_group_ids
                    .iter()
                    .any(|rg_id| !rg_to_sm.contains_key(rg_id));
            Ok((path.clone(), rg_to_sm, has_unmapped_read_groups))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut sample_names = input_mappings
        .iter()
        .flat_map(|(_, rg_to_sm, _)| rg_to_sm.values().cloned())
        .collect::<Vec<_>>();
    sample_names.sort();
    sample_names.dedup();

    if sample_names.is_empty() {
        let sample_names = vec![DEFAULT_SAMPLE_NAME.to_string()];
        let sample_index = sample_index(&sample_names);
        let input_resolvers = paths
            .iter()
            .cloned()
            .map(|source| InputSampleResolver {
                source,
                fallback_sample_index: Some(0),
                rg_to_sample_index: HashMap::new(),
                require_read_group: false,
            })
            .collect();
        return Ok(SampleResolution {
            sample_names,
            sample_index,
            input_resolvers,
        });
    }

    let sm_sample_count = sample_names.len();
    let fallback_input_indices = input_mappings
        .iter()
        .enumerate()
        .filter_map(|(input_index, (_, _, has_unmapped_read_groups))| {
            has_unmapped_read_groups.then_some(input_index)
        })
        .collect::<Vec<_>>();
    sample_names = unique_labels(
        sample_names.into_iter().chain(
            fallback_input_indices
                .iter()
                .map(|&input_index| file_label(&paths[input_index])),
        ),
    );
    let sample_index = sample_index(&sample_names);
    let mut fallback_sample_indices = vec![None; input_mappings.len()];
    for (fallback_offset, input_index) in fallback_input_indices.into_iter().enumerate() {
        fallback_sample_indices[input_index] = Some(sm_sample_count + fallback_offset);
    }
    let input_resolvers = input_mappings
        .into_iter()
        .zip(fallback_sample_indices)
        .map(
            |((source, rg_to_sm, _), fallback_sample_index)| InputSampleResolver {
                source,
                fallback_sample_index,
                rg_to_sample_index: rg_to_sm
                    .into_iter()
                    .map(|(rg, sm)| (rg, sample_index[&sm]))
                    .collect(),
                require_read_group: false,
            },
        )
        .collect();

    Ok(SampleResolution {
        sample_names,
        sample_index,
        input_resolvers,
    })
}

fn collect_file_resolution(paths: &[PathBuf]) -> Result<SampleResolution> {
    let sample_names = unique_file_labels(paths);
    let sample_index = sample_index(&sample_names);
    let input_resolvers = paths
        .iter()
        .cloned()
        .enumerate()
        .map(|(sample_index, source)| InputSampleResolver {
            source,
            fallback_sample_index: Some(sample_index),
            rg_to_sample_index: HashMap::new(),
            require_read_group: false,
        })
        .collect();

    Ok(SampleResolution {
        sample_names,
        sample_index,
        input_resolvers,
    })
}

fn collect_rg_resolution(paths: &[PathBuf]) -> Result<SampleResolution> {
    let file_labels = unique_file_labels(paths);
    let mut input_read_groups = Vec::with_capacity(paths.len());
    let mut raw_sample_names = Vec::new();

    for (path, file_label) in paths.iter().zip(&file_labels) {
        let mut read_groups = read_group_ids_from_header(&read_header(path)?)?;
        read_groups.sort();
        if read_groups.is_empty() {
            bail!(
                "--split-by rg requires at least one @RG header entry in {}",
                path.display()
            );
        }
        raw_sample_names.extend(
            read_groups
                .iter()
                .map(|rg| format!("{file_label}__{}", sanitize_vcf_sample_name(rg))),
        );
        input_read_groups.push(read_groups);
    }

    let sample_names = unique_labels(raw_sample_names);
    let sample_index = sample_index(&sample_names);
    let mut labels = sample_names.iter();
    let input_resolvers = paths
        .iter()
        .cloned()
        .zip(input_read_groups)
        .map(|(source, read_groups)| {
            let rg_to_sample_index = read_groups
                .into_iter()
                .map(|rg| {
                    let sample_name = labels
                        .next()
                        .expect("one sample name is generated for each read group");
                    (rg, sample_index[sample_name])
                })
                .collect();
            InputSampleResolver {
                source,
                fallback_sample_index: None,
                rg_to_sample_index,
                require_read_group: true,
            }
        })
        .collect();

    Ok(SampleResolution {
        sample_names,
        sample_index,
        input_resolvers,
    })
}

fn read_header(path: &Path) -> Result<sam::Header> {
    let file =
        File::open(path).with_context(|| format!("failed to open input BAM {}", path.display()))?;
    let mut reader = Reader::new(file);
    reader
        .read_header()
        .with_context(|| format!("failed to read header for {}", path.display()))
}

fn sample_index(sample_names: &[String]) -> HashMap<String, usize> {
    sample_names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.clone(), index))
        .collect()
}

fn read_group_ids_from_header(header: &sam::Header) -> Result<Vec<String>> {
    header
        .read_groups()
        .keys()
        .map(|rg_id| {
            std::str::from_utf8(rg_id.as_ref())
                .context("invalid RG id")
                .map(str::to_owned)
        })
        .collect()
}

fn unique_file_labels(paths: &[PathBuf]) -> Vec<String> {
    unique_labels(paths.iter().map(|path| file_label(path)))
}

fn file_label(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(sanitize_vcf_sample_name)
        .unwrap_or_else(|| "input".to_string())
}

fn sanitize_vcf_sample_name(value: &str) -> String {
    let sanitized = value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-') {
                character
            } else {
                '_'
            }
        })
        .collect::<String>();
    if sanitized.is_empty() {
        "input".to_string()
    } else {
        sanitized
    }
}

fn unique_labels(labels: impl IntoIterator<Item = String>) -> Vec<String> {
    let mut used = HashSet::new();
    let mut next_suffix = HashMap::new();
    let mut unique = Vec::new();

    for label in labels {
        if used.insert(label.clone()) {
            next_suffix.insert(label.clone(), 2usize);
            unique.push(label);
            continue;
        }

        let suffix = next_suffix.entry(label.clone()).or_insert(2);
        loop {
            let candidate = format!("{label}_{suffix}");
            *suffix += 1;
            if used.insert(candidate.clone()) {
                unique.push(candidate);
                break;
            }
        }
    }

    unique
}

pub(crate) fn read_rg_map(path: &Path) -> Result<Vec<(String, String)>> {
    let file =
        File::open(path).with_context(|| format!("failed to open rg_map {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    let mut saw_header = false;
    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("failed to read rg_map line")?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if !saw_header {
            if parts.len() != 2 || parts[0] != "RG" || parts[1] != "SM" {
                bail!("rg-map header must contain exactly RG and SM columns");
            }
            saw_header = true;
            continue;
        }
        if parts.len() != 2 {
            bail!(
                "rg-map row {} must contain exactly RG and SM columns",
                line_number + 1
            );
        }
        rows.push((parts[0].to_string(), parts[1].to_string()));
    }
    if !saw_header {
        bail!("rg-map header must contain exactly RG and SM columns");
    }
    Ok(rows)
}

fn build_rg_to_sm_from_map(
    rows: &[(String, String)],
) -> Result<(Vec<String>, HashMap<String, String>)> {
    let mut rg_to_sm = HashMap::new();
    let mut samples = Vec::new();

    for (rg, sm) in rows {
        insert_rg_sample(&mut rg_to_sm, rg.clone(), sm.clone(), "rg-map")?;
        if !samples.contains(sm) {
            samples.push(sm.clone());
        }
    }

    if rg_to_sm.is_empty() {
        bail!("rg-map did not contain any RG-to-sample mappings");
    }

    Ok((samples, rg_to_sm))
}

fn insert_rg_sample(
    rg_to_sm: &mut HashMap<String, String>,
    rg: String,
    sm: String,
    source: &str,
) -> Result<()> {
    if let Some(existing) = rg_to_sm.get(&rg) {
        if existing != &sm {
            bail!(
                "conflicting sample mapping for RG {}: {} vs {} ({})",
                rg,
                existing,
                sm,
                source
            );
        }
        return Ok(());
    }
    rg_to_sm.insert(rg, sm);
    Ok(())
}

fn rg_to_sm_from_header(header: &sam::Header) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for (rg_id, read_group) in header.read_groups() {
        let rg_id = std::str::from_utf8(rg_id.as_ref()).context("invalid RG id")?;
        let sm = read_group
            .other_fields()
            .get(&sam::header::record::value::map::read_group::tag::SAMPLE)
            .map(|value| std::str::from_utf8(value.as_ref()).context("invalid SM value"))
            .transpose()?;
        if let Some(sm) = sm {
            map.insert(rg_id.to_string(), sm.to_string());
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::{
        DEFAULT_SAMPLE_NAME, build_rg_to_sm_from_map, read_rg_map, unique_file_labels,
        unique_labels,
    };
    use anyhow::Result;
    use std::{io::Write, path::PathBuf};
    use tempfile::NamedTempFile;

    #[test]
    fn read_rg_map_parses_tab_separated_lines() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "RG\tSM")?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "rg2\tsample_b")?;

        let rows = read_rg_map(f.path())?;
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0], ("rg1".to_string(), "sample_a".to_string()));
        assert_eq!(rows[1], ("rg2".to_string(), "sample_b".to_string()));
        Ok(())
    }

    #[test]
    fn read_rg_map_skips_comments_and_blank_lines() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "# header comment")?;
        writeln!(f)?;
        writeln!(f, "RG\tSM")?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "rg2\tsample_b")?;

        let rows = read_rg_map(f.path())?;
        assert_eq!(rows.len(), 2);
        Ok(())
    }

    #[test]
    fn read_rg_map_errors_on_malformed_rows() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "RG\tSM")?;
        writeln!(f, "rg1\tsample_a")?;
        writeln!(f, "orphan_rg")?; // malformed row
        writeln!(f, "rg2\tsample_b")?;

        let err = read_rg_map(f.path()).unwrap_err();
        assert!(err.to_string().contains("exactly RG and SM columns"));
        Ok(())
    }

    #[test]
    fn read_rg_map_errors_without_rg_sm_header() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "rg1\tsample_a")?;

        let err = read_rg_map(f.path()).unwrap_err();
        assert!(
            err.to_string()
                .contains("header must contain exactly RG and SM")
        );
        Ok(())
    }

    #[test]
    fn read_rg_map_errors_on_empty_file() {
        let f = NamedTempFile::new().unwrap();

        let err = read_rg_map(f.path()).unwrap_err();
        assert!(
            err.to_string()
                .contains("header must contain exactly RG and SM")
        );
    }

    #[test]
    fn build_rg_to_sm_preserves_insertion_order_for_samples() -> Result<()> {
        let rows = vec![
            ("rg1".to_string(), "sample_a".to_string()),
            ("rg2".to_string(), "sample_b".to_string()),
            ("rg3".to_string(), "sample_b".to_string()), // second RG for same sample
        ];
        let (samples, rg_to_sm) = build_rg_to_sm_from_map(&rows)?;
        assert_eq!(samples, vec!["sample_a", "sample_b"]);
        assert_eq!(rg_to_sm["rg1"], "sample_a");
        assert_eq!(rg_to_sm["rg3"], "sample_b");
        Ok(())
    }

    #[test]
    fn build_rg_to_sm_errors_on_conflicting_mapping() {
        let rows = vec![
            ("rg1".to_string(), "sample_a".to_string()),
            ("rg1".to_string(), "sample_b".to_string()), // same RG, different sample
        ];
        let err = build_rg_to_sm_from_map(&rows).unwrap_err();
        assert!(err.to_string().contains("conflicting sample mapping"));
    }

    #[test]
    fn build_rg_to_sm_errors_on_empty_input() {
        let rows: Vec<(String, String)> = vec![];
        let err = build_rg_to_sm_from_map(&rows).unwrap_err();
        assert!(err.to_string().contains("did not contain any"));
    }

    #[test]
    fn default_sample_name_is_vcf_safe() {
        assert_eq!(DEFAULT_SAMPLE_NAME, "SAMPLE");
    }

    #[test]
    fn unique_labels_adds_stable_suffixes() {
        let labels = unique_labels(
            ["sample", "sample_2", "sample"]
                .into_iter()
                .map(str::to_owned),
        );
        assert_eq!(labels, ["sample", "sample_2", "sample_3"]);
    }

    #[test]
    fn file_labels_are_vcf_safe_and_collision_free() {
        let labels =
            unique_file_labels(&[PathBuf::from("run one.bam"), PathBuf::from("run one.cram")]);
        assert_eq!(labels, ["run_one", "run_one_2"]);
    }
}
