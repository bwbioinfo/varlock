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
- `annotate` - VCF annotation from VCF-like databases.
- `filter` - INFO-based VCF filtering.
- `intersect` - exact-key VCF intersection and difference.

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
cargo run -- filter --help
cargo run -- intersect --help
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

## Paired Calling

Use `--pair tumor=SAMPLE,normal=SAMPLE` to filter calls through one tumor/normal
pair after the global ALT has been selected. The tumor and normal names must
match sample names from BAM `RG/SM` tags or `--rg-map`.

```bash
varlock call-targets \
  --input tumor.bam \
  --input normal.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --pair tumor=Tumor,normal=Normal \
  --tumor-min-alt-count 3 \
  --tumor-min-alt-fraction 0.05 \
  --normal-max-alt-count 0 \
  --normal-max-alt-fraction 0.0 \
  --normal-min-depth 10 \
  --output paired.calls.vcf.gz
```

Paired output keeps all samples in `FORMAT/GT:DP:AD` and adds paired INFO fields
to records that pass:

- `PAIR` - tumor and normal sample names.
- `SOMATIC` - flag indicating the selected ALT passed paired filters.
- `TUMOR_AF`, `NORMAL_AF` - selected ALT allele fractions.
- `TUMOR_ALT_COUNT`, `NORMAL_ALT_COUNT` - selected ALT counts.
- `TUMOR_DP`, `NORMAL_DP` - tumor and normal depths.

GPU calling accepts the same paired options and uses the same output
classification:

```bash
varlock call-targets-gpu \
  --input tumor.bam \
  --input normal.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --pair tumor=Tumor,normal=Normal \
  --tumor-min-alt-count 3 \
  --normal-max-alt-count 0 \
  --output paired.gpu.calls.vcf.gz
```

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

`annotate` adds INFO annotations to an input VCF using `CHROM, POS, REF, ALT`
matches from one or more VCF-like databases. By default, matching is exact. Use
`--reference` to normalize simple indel representations before matching; this
preserves the original input records while comparing normalized keys. If no
reference is supplied and a database has an adjacent `.tbi` or `.csi` index,
varlock queries it by input variant position; otherwise it loads that database
into memory. Normalized matching currently loads databases into memory so their
keys can be normalized. Multi-ALT database records are split by ALT for matching,
and mapped comma-valued fields with one value per ALT are emitted in input ALT
order.

```bash
varlock annotate \
  --input calls.vcf.gz \
  --database gnomad=gnomad.sites.vcf.gz \
  --annotation gnomad:AF=gnomAD_AF,AC=gnomAD_AC \
  --reference hg38/hg38.fa \
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

## VCF Filtering

`filter` streams a VCF and keeps records that pass simple INFO predicates. This
filtering mode is intended for post-annotation frequency, marker-field, and
sample-aware FORMAT filters.

```bash
varlock filter \
  --input calls.annotated.vcf.gz \
  --require-info gnomAD_AF \
  --max-info gnomAD_AF=0.01 \
  --output calls.rare.vcf.gz
```

Records can also be dropped when an INFO field is present:

```bash
varlock filter \
  --input calls.annotated.vcf.gz \
  --exclude-info common_AF \
  --output calls.no_common.vcf.gz
```

For comma-valued INFO fields, `--max-info FIELD=VALUE` requires every numeric
non-missing value in the field to be at or below the threshold.

For more complex INFO predicates, use `--expr`. Repeated expressions are
combined with logical AND:

```bash
varlock filter \
  --input calls.annotated.vcf.gz \
  --expr "gnomAD_AF < 0.01 && DP >= 20" \
  --expr "CLNSIG != 'Benign' && missing(COMMON)" \
  --output calls.expr_filtered.vcf.gz
```

The first-pass expression language supports INFO field names, numeric and quoted
string literals, `==`, `!=`, `<`, `<=`, `>`, `>=`, `&&`, `||`, `!`,
parentheses, and `missing(FIELD)`. A bare INFO field name is true when that
field is present, so `SOMATIC && !COMMON` is valid.

Sample-aware filters inspect `FORMAT` values for named samples:

```bash
varlock filter \
  --input paired.calls.vcf.gz \
  --sample-has-alt Tumor \
  --sample-min-dp Tumor=10 \
  --sample-gt Normal=0/0 \
  --output paired.somatic_like.vcf.gz
```

Groups can be defined once and reused by group predicates:

```bash
varlock filter \
  --input cohort.calls.vcf.gz \
  --sample-group affected=TumorA,TumorB \
  --sample-group controls=NormalA,NormalB \
  --group-any-has-alt affected \
  --group-all-min-dp controls=20 \
  --output cohort.group_filtered.vcf.gz
```

Available first-pass sample predicates are:

- `--sample-has-alt SAMPLE` - `FORMAT/GT` contains any non-reference allele.
- `--sample-gt SAMPLE=GT` - exact `FORMAT/GT` match.
- `--sample-min-dp SAMPLE=DP` - `FORMAT/DP` is at least `DP`.
- `--group-any-has-alt GROUP` - any group member has an alternate allele.
- `--group-any-gt GROUP=GT` - any group member has exact `FORMAT/GT`.
- `--group-all-min-dp GROUP=DP` - every group member has `FORMAT/DP` at least
  `DP`.

Filter output preserves all input samples and is bgzipped VCF with a CSI index
by default; use `--index-type tbi` to write a tabix index instead.

## VCF Intersection And Difference

`intersect` performs exact allele-key set operations between two VCFs, or
genotype-aware set operations between two sample groups in one multi-sample VCF.
Two-file mode matches variants by `(CHROM, POS, REF, ALT)` and treats record
presence as support. Multi-ALT records are split internally for matching; output
preserves the original record from the emitted side.

Shared variants emit records from the left VCF:

```bash
varlock intersect \
  --left tumor.vcf.gz \
  --right normal.vcf.gz \
  --mode shared \
  --output shared.vcf.gz
```

Difference modes emit records unique to one side:

```bash
varlock intersect \
  --left tumor.vcf.gz \
  --right normal.vcf.gz \
  --mode left-only \
  --output tumor_only.vcf.gz
```

```bash
varlock intersect \
  --left tumor.vcf.gz \
  --right normal.vcf.gz \
  --mode right-only \
  --output normal_only.vcf.gz
```

Output records include `INFO/VARLOCK_SET=both`, `left`, or `right` and are
written as bgzipped VCF with a CSI index by default.

One-file sample-group mode compares genotype support between explicit sample
groups. A sample supports a variant when `FORMAT/GT` contains any non-reference
allele. Missing genotypes such as `./.` or `.|.`, missing alleles, and records
without `FORMAT/GT` count as absence. Output preserves all input samples and adds
`INFO/VARLOCK_LEFT_SUPPORT` and `INFO/VARLOCK_RIGHT_SUPPORT`.

```bash
varlock intersect \
  --input cohort.vcf.gz \
  --left-samples tumor_a,tumor_b \
  --right-samples normal_a,normal_b \
  --mode left-only \
  --output tumor_group_only.vcf.gz
```

## Roadmap

Planned features are tracked with `bd` issues.

### Variant Intersection And Difference

Implemented:

- two-file exact-key `shared`, `left-only`, and `right-only`
- one-file multi-sample group `shared`, `left-only`, and `right-only`

Planned examples:

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

Planned implementation phases:

1. Add named multi-file sets with repeated `--set NAME=vcf.gz:s1,s2` and a
   manifest form for `all-shared`, `any-shared`, and `set-diff A-B`.
2. Add matching-mode flags: `--reference` for normalized allele keys,
   `--site-only`, and explicit `--genotype-aware`.

### Annotation Follow-Ups

The annotation command supports exact and reference-normalized VCF-like database
matching, multi-ALT records, and indexed lookup for `.tbi` and `.csi` databases.
Future work should add TSV/BED-style variant or interval database adapters.

Design decisions to settle before implementation:

- required database indexing and supported database formats
- exact allele matching vs site-only annotation
- optional normalization against a reference
- INFO/header naming and provenance lines
- no-hit representation
- multi-database conflict handling
