# varlock

`varlock` is a Rust CLI for simple pileup-based SNV calling from BAM alignments.
It writes bgzipped VCF output and an index. The current implementation focuses on
calling A/C/G/T substitutions from one or more BAMs, either across target regions
or across all covered positions.

## Current Commands

```bash
varlock --help
```

Current subcommands:

- `call-targets` - CPU SNV calling from BAMs.
- `call-targets-gpu` - GPU-accelerated count aggregation, available when built
  with the `wgpu` feature.
- `annotate` - exact-key VCF annotation from VCF-like databases.

## Build

```bash
cargo build
cargo test
```

Build with GPU support:

```bash
cargo build --release --features wgpu
cargo test --features wgpu
```

During development:

```bash
cargo run -- call-targets --help
cargo run --features wgpu -- call-targets-gpu --help
cargo run -- annotate --help
```

## Inputs

`varlock` accepts BAM input paths through either repeated `--input` arguments or
a `--bamlist` file.

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output sample.calls.vcf.gz
```

Directories passed to `--input` are searched recursively for `.bam` files.

```bash
varlock call-targets \
  --input ./bams \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output cohort.calls.vcf.gz
```

Use a BAM list when the input set is easier to describe in a file:

```bash
varlock call-targets \
  --bamlist bams.txt \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output cohort.calls.vcf.gz
```

`bams.txt` contains one BAM path per line. Blank lines and lines beginning with
`#` are skipped.

## Reference FASTA

The reference may be plain FASTA or bgzipped FASTA.

```bash
--reference hg38/hg38.fa
```

If the reference is missing the required FASTA index, `varlock` prepares a sorted
and indexed reference beside the input reference.

## Targets And No-Target Mode

Use `--targets` with a BED file for target-region calling. BED intervals are
0-based, half-open.

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets panel.bed
```

Omit `--targets` to call all covered positions:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --output sample.covered.vcf.gz
```

## Samples And Read Groups

By default, `varlock` uses BAM read group `SM` tags to determine samples. If a
BAM has no read group/sample mapping, it is treated as one default sample.

Provide an explicit read-group to sample map with `--rg-map`:

```bash
varlock call-targets \
  --input tumor.bam \
  --input normal.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --rg-map rg_to_sample.tsv \
  --output paired.calls.vcf.gz
```

`rg_to_sample.tsv` format:

```text
RG001    Tumor
RG002    Normal
```

## Calling Filters

Common filters:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --min-mapq 20 \
  --min-baseq 20 \
  --min-alt-count 3 \
  --min-alt-fraction 0.05 \
  --max-depth 100000
```

`--min-alt-count` and `--min-alt-fraction` are evaluated across all samples at a
site. The selected ALT is the highest-count non-reference A/C/G/T base after
summing counts across samples. Per-sample genotype fields are still written for
the selected ALT.

## Output

Output is a bgzipped VCF:

```bash
--output calls.vcf.gz
```

The default index type is CSI. TBI is also available:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --index-type tbi
```

VCF records include:

- `INFO/DP` - total depth across samples.
- `FORMAT/GT:DP:AD` - genotype, sample depth, and reference/alternate allele
  depths for each sample.

## GPU Calling

Build with `wgpu`:

```bash
cargo build --release --features wgpu
```

Run GPU calling:

```bash
varlock call-targets-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output sample.gpu.vcf.gz
```

If no compatible GPU is found, `call-targets-gpu` falls back to CPU calling
unless `--require-gpu` is set:

```bash
varlock call-targets-gpu \
  --require-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

No-target GPU calling is supported:

```bash
varlock call-targets-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --output sample.gpu.covered.vcf.gz
```

No-target GPU mode flushes covered observations in batches, merges raw counts,
and applies `--max-depth` once before output so results do not depend on flush
threshold. The flush threshold can be tuned:

```bash
varlock call-targets-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --obs-flush-threshold 500000
```

## GPU Selection And Multi-GPU Static Targets

List compatible adapters:

```bash
varlock call-targets-gpu --gpu-list
```

Select one GPU by index:

```bash
varlock -v call-targets-gpu \
  --gpu-index 0 \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

Select one GPU by name substring:

```bash
varlock call-targets-gpu \
  --gpu-name H100 \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

Use multiple GPUs for static target calling by repeating `--gpu-index`:

```bash
varlock -v call-targets-gpu \
  --gpu-index 0 \
  --gpu-index 1 \
  --gpu-index 2 \
  --gpu-index 3 \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output sample.multi-gpu.vcf.gz
```

Current multi-GPU execution applies to the static target path. Include
`--targets` to shard target chunks across selected GPU adapters. No-target GPU
calling is correctness-fixed but is not currently sharded across multiple GPUs.

GPU tuning options:

```bash
varlock call-targets-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --matrix-budget-mib 256 \
  --max-obs-upload 1000000 \
  --obs-flush-threshold 500000
```

## Logging

Use global verbosity flags before the subcommand:

```bash
varlock -v call-targets-gpu --input sample.bam --reference hg38/hg38.fa --targets targets.bed
varlock -vv call-targets --input sample.bam --reference hg38/hg38.fa --targets targets.bed
```

Mirror stderr logs to a file:

```bash
varlock --log-file run.log call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

## Variant Annotation

`annotate` adds INFO annotations to an input VCF using exact
`CHROM, POS, REF, ALT` matches from one or more VCF-like databases. If a
database has an adjacent `.tbi` or `.csi` index, varlock queries it by input
variant position; otherwise it loads that database into memory. The current
annotation mode supports single-ALT records. Multi-ALT database and input
records are left unannotated.

```bash
varlock annotate \
  --input calls.vcf.gz \
  --database gnomad=gnomad.sites.vcf.gz \
  --annotation gnomad:AF=gnomAD_AF,AC=gnomAD_AC \
  --index-type csi \
  --output calls.annotated.vcf.gz
```

Multiple databases can be supplied by repeating `--database` and `--annotation`:

```bash
varlock annotate \
  --input calls.vcf.gz \
  --database common=dbs/common.vcf.gz \
  --database cohort=dbs/cohort_freqs.vcf.gz \
  --annotation common:AF=common_AF \
  --annotation cohort:AF=cohort_AF \
  --output calls.annotated.vcf.gz
```

Annotation output is bgzipped VCF. Header lines are added for each destination
INFO field and for each annotation database. Output indexes are written by
default as CSI (`calls.annotated.vcf.gz.csi`); use `--index-type tbi` to write a
tabix index instead.

## Roadmap

Planned features are tracked with `bd` issues.

### Variant Intersection And Difference

Planned command family for set operations over VCFs:

- two single-sample VCFs
- one multi-sample VCF with explicit sample groups
- multiple multi-sample VCFs with explicit file/sample sets

Proposed examples:

```bash
varlock intersect \
  --left tumor.vcf.gz \
  --right normal.vcf.gz \
  --mode shared \
  --output shared.vcf.gz
```

```bash
varlock intersect \
  --input cohort.vcf.gz \
  --left-samples a,b \
  --right-samples c,d \
  --mode left-only \
  --output ab_not_cd.vcf.gz
```

```bash
varlock intersect \
  --set A=file1.vcf.gz:a,b \
  --set B=file2.vcf.gz:c,d \
  --mode set-diff A-B \
  --output A_minus_B.vcf.gz
```

Design decisions to settle before implementation:

- variant identity: exact allele key `(CHROM, POS, REF, ALT)` vs site-only
- optional left-normalization against a reference
- genotype-aware matching vs presence/absence matching
- whether missing genotypes count as absence
- projected VCF output vs tabular summaries

### Annotation Follow-Ups

The initial annotation command supports exact single-ALT VCF-like databases and
indexed lookup for `.tbi` and `.csi` databases. Future work should add
normalization, multi-ALT handling, and TSV/BED-style variant or interval database
adapters.

Design decisions to settle before implementation:

- required database indexing and supported database formats
- exact allele matching vs site-only annotation
- optional normalization against a reference
- INFO/header naming and provenance lines
- no-hit representation
- multi-database conflict handling
