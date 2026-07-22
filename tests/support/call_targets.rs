use std::{
    fs::File,
    io::{BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::{Cigar, Data, QualityScores, RecordBuf, Sequence, data::field::Value},
    },
    header::record::value::{
        Map,
        map::{
            Header, ReadGroup, ReferenceSequence,
            header::{sort_order::COORDINATE, tag::SORT_ORDER},
            read_group::tag::SAMPLE,
        },
    },
};
use noodles_vcf as vcf;

const FASTA_LINE_BASES: usize = 60;

#[derive(Clone, Debug)]
pub struct ReadSpec {
    name: String,
    reference_sequence_id: usize,
    alignment_start: usize,
    sequence: Vec<u8>,
    mapping_quality: u8,
    base_quality: u8,
    read_group: Option<String>,
}

impl ReadSpec {
    pub fn single_base(name: impl Into<String>, alignment_start: usize, base: u8) -> Self {
        Self {
            name: name.into(),
            reference_sequence_id: 0,
            alignment_start,
            sequence: vec![base],
            mapping_quality: 60,
            base_quality: 40,
            read_group: None,
        }
    }

    pub fn with_read_group(mut self, read_group: impl Into<String>) -> Self {
        self.read_group = Some(read_group.into());
        self
    }

    #[allow(dead_code)]
    pub fn on_reference(mut self, reference_sequence_id: usize) -> Self {
        self.reference_sequence_id = reference_sequence_id;
        self
    }
}

#[derive(Clone, Debug)]
struct ReferenceSpec {
    name: String,
    length: usize,
}

#[derive(Clone, Debug)]
struct ReadGroupSpec {
    id: String,
    sample: Option<String>,
}

#[derive(Debug)]
pub struct BamFixture {
    pub path: PathBuf,
    pub bai_path: PathBuf,
}

#[derive(Debug)]
pub struct BamFixtureBuilder {
    references: Vec<ReferenceSpec>,
    read_groups: Vec<ReadGroupSpec>,
    records: Vec<ReadSpec>,
}

impl BamFixtureBuilder {
    pub fn new(reference_name: impl Into<String>, reference_length: usize) -> Self {
        Self {
            references: vec![ReferenceSpec {
                name: reference_name.into(),
                length: reference_length,
            }],
            read_groups: Vec::new(),
            records: Vec::new(),
        }
    }

    #[allow(dead_code)]
    pub fn add_reference(&mut self, name: impl Into<String>, length: usize) {
        self.references.push(ReferenceSpec {
            name: name.into(),
            length,
        });
    }

    pub fn add_read_group(&mut self, id: impl Into<String>, sample: Option<&str>) {
        self.read_groups.push(ReadGroupSpec {
            id: id.into(),
            sample: sample.map(str::to_owned),
        });
    }

    pub fn add_read(&mut self, read: ReadSpec) {
        self.records.push(read);
    }

    pub fn add_observations(
        &mut self,
        alignment_start: usize,
        base: u8,
        count: usize,
        read_group: Option<&str>,
    ) {
        for _ in 0..count {
            let mut read = ReadSpec::single_base(
                format!("read-{:06}", self.records.len()),
                alignment_start,
                base,
            );
            if let Some(read_group) = read_group {
                read = read.with_read_group(read_group);
            }
            self.add_read(read);
        }
    }

    pub fn write(self, path: impl AsRef<Path>) -> Result<BamFixture> {
        let path = path.as_ref();
        let header = self.build_header()?;
        let file = File::create(path)
            .with_context(|| format!("failed to create BAM fixture {}", path.display()))?;
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header)?;

        for read in self.records {
            let record = build_record(read)?;
            writer.write_alignment_record(&header, &record)?;
        }
        writer.try_finish()?;
        drop(writer);

        let index = bam::fs::index(path)
            .with_context(|| format!("failed to index BAM fixture {}", path.display()))?;
        let bai_path = append_extension(path, "bai");
        bam::bai::fs::write(&bai_path, &index)
            .with_context(|| format!("failed to write BAI fixture {}", bai_path.display()))?;

        Ok(BamFixture {
            path: path.to_path_buf(),
            bai_path,
        })
    }

    fn build_header(&self) -> Result<sam::Header> {
        let header_map = Map::<Header>::builder()
            .insert(SORT_ORDER, COORDINATE)
            .build()?;
        let mut builder = sam::Header::builder().set_header(header_map);

        for reference in &self.references {
            let length = NonZero::new(reference.length)
                .with_context(|| format!("reference {} has zero length", reference.name))?;
            builder = builder.add_reference_sequence(
                reference.name.as_str(),
                Map::<ReferenceSequence>::new(length),
            );
        }

        for read_group in &self.read_groups {
            let map = match read_group.sample.as_deref() {
                Some(sample) => Map::<ReadGroup>::builder().insert(SAMPLE, sample).build()?,
                None => Map::<ReadGroup>::default(),
            };
            builder = builder.add_read_group(read_group.id.as_str(), map);
        }

        Ok(builder.build())
    }
}

fn build_record(read: ReadSpec) -> Result<RecordBuf> {
    if read.sequence.is_empty() {
        bail!("BAM fixture read {} has an empty sequence", read.name);
    }

    let alignment_start = Position::try_from(read.alignment_start)
        .with_context(|| format!("invalid alignment start for fixture read {}", read.name))?;
    let mapping_quality = MappingQuality::new(read.mapping_quality)
        .with_context(|| format!("invalid mapping quality for fixture read {}", read.name))?;
    let cigar: Cigar = [Op::new(Kind::Match, read.sequence.len())]
        .into_iter()
        .collect();
    let quality_scores = QualityScores::from(vec![read.base_quality; read.sequence.len()]);
    let mut data = Data::default();
    if let Some(read_group) = read.read_group {
        data.insert(Tag::READ_GROUP, Value::from(read_group));
    }

    Ok(RecordBuf::builder()
        .set_name(read.name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(read.reference_sequence_id)
        .set_alignment_start(alignment_start)
        .set_mapping_quality(mapping_quality)
        .set_cigar(cigar)
        .set_sequence(Sequence::from(read.sequence))
        .set_quality_scores(quality_scores)
        .set_data(data)
        .build())
}

#[derive(Debug)]
pub struct FastaFixture {
    pub path: PathBuf,
    pub fai_path: PathBuf,
}

pub fn write_fasta(path: impl AsRef<Path>, records: &[(&str, &[u8])]) -> Result<FastaFixture> {
    if records.is_empty() {
        bail!("FASTA fixture requires at least one reference sequence");
    }

    let path = path.as_ref();
    let fai_path = append_extension(path, "fai");
    let mut fasta = BufWriter::new(
        File::create(path)
            .with_context(|| format!("failed to create FASTA fixture {}", path.display()))?,
    );
    let mut fai = BufWriter::new(
        File::create(&fai_path)
            .with_context(|| format!("failed to create FAI fixture {}", fai_path.display()))?,
    );
    let mut offset = 0u64;

    for &(name, sequence) in records {
        if sequence.is_empty() {
            bail!("FASTA fixture reference {name} has zero length");
        }

        let header = format!(">{name}\n");
        fasta.write_all(header.as_bytes())?;
        offset += header.len() as u64;
        let sequence_offset = offset;

        for chunk in sequence.chunks(FASTA_LINE_BASES) {
            fasta.write_all(chunk)?;
            fasta.write_all(b"\n")?;
            offset += (chunk.len() + 1) as u64;
        }

        let line_bases = sequence.len().min(FASTA_LINE_BASES);
        writeln!(
            fai,
            "{name}\t{}\t{sequence_offset}\t{line_bases}\t{}",
            sequence.len(),
            line_bases + 1
        )?;
    }

    fasta.flush()?;
    fai.flush()?;

    Ok(FastaFixture {
        path: path.to_path_buf(),
        fai_path,
    })
}

pub fn write_bed(path: impl AsRef<Path>, intervals: &[(&str, u64, u64)]) -> Result<PathBuf> {
    let path = path.as_ref();
    let mut writer = BufWriter::new(
        File::create(path)
            .with_context(|| format!("failed to create BED fixture {}", path.display()))?,
    );
    for &(reference, start, end) in intervals {
        writeln!(writer, "{reference}\t{start}\t{end}")?;
    }
    writer.flush()?;
    Ok(path.to_path_buf())
}

pub struct ParsedVcf {
    pub header: vcf::Header,
    pub records: Vec<vcf::variant::RecordBuf>,
}

pub fn read_vcf(path: impl AsRef<Path>) -> Result<ParsedVcf> {
    let path = path.as_ref();
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("failed to open VCF fixture output {}", path.display()))?;
    let header = reader.read_header()?;
    let records = reader
        .records()
        .map(|result| {
            let record = result?;
            vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
        })
        .collect::<std::io::Result<Vec<_>>>()?;

    Ok(ParsedVcf { header, records })
}

pub fn append_extension(path: impl AsRef<Path>, extension: &str) -> PathBuf {
    let mut value = path.as_ref().as_os_str().to_os_string();
    value.push(".");
    value.push(extension);
    PathBuf::from(value)
}
