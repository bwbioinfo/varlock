# varlock

`varlock` is a Rust CLI for pileup-based SNV calling, VCF annotation, filtering,
and variant set comparison. The calling workflow reads one or more BAM files,
counts A/C/G/T support at target sites or all covered positions, and writes a
bgzipped VCF plus an index.

The caller is designed for targeted panels and cohort-style workflows. Sample
identity can come from BAM read group `SM` tags, input BAM files, input-local
read groups, or an explicit read-group map. With the optional `wgpu` feature,
`call-targets` uses GPU aggregation by default and falls back to the CPU path
when no compatible adapter is available.

## Quick Start

Build and run the test suite:

```bash
cargo build --release
cargo test
```

Call variants in a targeted BAM:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets panel.bed \
  --output sample.calls.vcf.gz
```

Build with GPU support and use the same command:

```bash
cargo build --release --features wgpu
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets panel.bed \
  --output sample.calls.vcf.gz
```

Force CPU execution in a GPU-enabled build:

```bash
varlock call-targets --cpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets panel.bed
```

## Current Commands

```bash
varlock --help
```

Current subcommands:

- `call-targets` - SNV calling from BAMs. When built with the `wgpu` feature,
  it uses GPU acceleration by default and falls back to CPU if no compatible GPU
  is found. Pass `--cpu` to force CPU calling.
- `call-targets-gpu` - compatibility alias for GPU-accelerated count
  aggregation, available when built with the `wgpu` feature.
- `annotate` - VCF annotation from VCF-like databases.
- `filter` - INFO-based VCF filtering.
- `intersect` - exact-key VCF intersection and difference.

Most commands write bgzipped VCF output and create a CSI index by default. Use
`--index-type tbi` when a tabix index is required.

## Requirements

- Rust toolchain with Cargo.
- Input BAM files with coordinate-sorted alignments and readable headers.
- Reference FASTA matching the BAM reference names.
- Optional BED targets for targeted calling.
- Optional GPU support through the `wgpu` Cargo feature.

## Build

```bash
cargo build
cargo test
```

Build optimized binaries:

```bash
cargo build --release
```

Build and test with GPU support:

```bash
cargo build --release --features wgpu
cargo test --features wgpu
```

During development:

```bash
cargo run -- call-targets --help
cargo run --features wgpu -- call-targets --help
cargo run -- annotate --help
cargo run -- filter --help
cargo run -- intersect --help
```

## Caller Overview

`call-targets` performs simple SNV calling for A/C/G/T substitutions. It does
not currently call indels, symbolic variants, or structural variants. The caller
chooses the highest-count non-reference base at each retained site after
summing support across samples, then writes per-sample `GT:DP:AD` fields for the
selected ALT.

Core calling inputs:

| Input | Purpose |
|---|---|
| `--input PATH` | BAM file or directory of BAM files; repeat for multiple paths. |
| `--bamlist FILE` | Text file with one BAM path per line. |
| `--reference FASTA` | Reference FASTA used for sequence names and REF bases. |
| `--targets BED` | Optional target intervals; omit to call all covered positions. |
| `--split-by MODE` | VCF sample identity policy: `sm` (default), `file`, or `rg`. |
| `--rg-map FILE` | Optional read-group to sample TSV with `RG` and `SM` headers. |

Core calling thresholds:

| Option | Meaning |
|---|---|
| `--min-mapq` | Minimum read mapping quality. |
| `--min-baseq` | Minimum base quality. |
| `--min-alt-count` | Minimum selected ALT count across samples. |
| `--min-alt-fraction` | Minimum selected ALT fraction across samples. |
| `--max-depth` | Per-site depth cap used while accumulating counts. |

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

Reference sequence names must match the BAM header sequence names. If the
reference is prepared automatically, the prepared FASTA is reused on later runs
when the expected index files are already present.

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

`--split-by sm` is the default and preserves the existing SM-based behavior:
read groups with the same BAM-header `SM` value share a VCF sample column. If
every input lacks `SM`, all reads use one `SAMPLE` column. Reused RG IDs are
resolved per input BAM, so the same RG ID can map to different `SM` values in
different BAMs. If only some inputs lack `SM`, reads without an SM mapping use
that input's VCF-safe filename stem rather than being dropped.

Use `--split-by file` for one VCF sample column per BAM, regardless of header
tags. This is the recommended choice when each input BAM is one biological
sample but its `@RG` records omit `SM`:

```bash
varlock call-targets \
  --input tumor.bam \
  --input normal.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --split-by file \
  --output paired.calls.vcf.gz
```

The resulting sample labels are VCF-safe filename stems, such as `tumor` and
`normal`; duplicate stems receive deterministic suffixes (`sample`,
`sample_2`).

Use `--split-by rg` when lanes or libraries need independent columns. Each
column is named `file-stem__rg-id`, so repeated RG IDs from different BAMs stay
distinct. Every retained read must have an RG tag declared by that input BAM;
missing or unknown tags stop the command before output is written.

```bash
varlock call-targets \
  --input cohort_a.bam \
  --input cohort_b.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --split-by rg \
  --output lanes.calls.vcf.gz
```

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
RG       SM
RG001    Tumor
RG002    Normal
```

Use `--rg-map` with `--split-by sm` when BAM headers have incorrect sample
names or need to be overridden for a particular analysis. The TSV must have
exactly `RG` and `SM` headers and two fields on every data row. `--rg-map` is
not accepted with `--split-by file` or `--split-by rg`, because a two-column
map cannot disambiguate duplicate RG IDs across inputs.

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

For tumor/normal work, use paired calling thresholds in addition to the global
thresholds so a globally selected ALT must also satisfy tumor and normal sample
constraints.

## Paired Calling

Use `--pair tumor=SAMPLE,normal=SAMPLE` to filter calls through one tumor/normal
pair after the global ALT has been selected. The tumor and normal names must
match resolved sample labels from `--split-by sm`, `--split-by file`, or
`--split-by rg`.

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

GPU calling accepts the same paired options through `call-targets` when built
with `wgpu` and uses the same output classification:

```bash
varlock call-targets \
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

For a selected ALT, `AD` is emitted as reference depth followed by selected ALT
depth. Sites that do not satisfy the global filters, paired filters, or basic
A/C/G/T constraints are not emitted.

Output paths:

- explicit `--output` paths are used as provided
- without `--output`, the output name is derived from the first input or BAM
  list path
- CSI indexes are written as `<output>.csi`
- TBI indexes are written as `<output>.tbi`

## GPU Calling

Build with `wgpu`:

```bash
cargo build --release --features wgpu
```

With `wgpu` enabled, `call-targets` tries GPU acceleration by default:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --output sample.gpu.vcf.gz
```

Force the CPU path with `--cpu`:

```bash
varlock call-targets \
  --cpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

If no compatible GPU is found, `call-targets` falls back to CPU calling unless
`--require-gpu` is set:

```bash
varlock call-targets \
  --require-gpu \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

No-target GPU calling is supported:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --output sample.gpu.covered.vcf.gz
```

No-target GPU mode flushes covered observations in batches, merges raw counts,
and applies `--max-depth` once before output so results do not depend on flush
threshold. The flush threshold can be tuned:

```bash
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --obs-flush-threshold 500000
```

GPU and CPU paths share the same sample collection, target loading, paired
calling, and VCF output code. The GPU path accelerates count aggregation; it
does not change the calling thresholds or output schema.

## GPU Selection And Multi-GPU Static Targets

List compatible adapters:

```bash
varlock call-targets --gpu-list
```

Select one GPU by index:

```bash
varlock -v call-targets \
  --gpu-index 0 \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

Select one GPU by name substring:

```bash
varlock call-targets \
  --gpu-name H100 \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

Use multiple GPUs for static target calling by repeating `--gpu-index`:

```bash
varlock -v call-targets \
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
varlock call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed \
  --matrix-budget-mib 256 \
  --max-obs-upload 1000000 \
  --obs-flush-threshold 500000
```

GPU options are available on `call-targets` only when the binary is built with
`--features wgpu`. The legacy `call-targets-gpu` subcommand remains available as
a compatibility alias in GPU-enabled builds, but new examples should prefer
`call-targets`.

## Logging

Use global verbosity flags before the subcommand:

```bash
varlock -v call-targets --input sample.bam --reference hg38/hg38.fa --targets targets.bed
varlock -vv call-targets --cpu --input sample.bam --reference hg38/hg38.fa --targets targets.bed
```

Mirror stderr logs to a file:

```bash
varlock --log-file run.log call-targets \
  --input sample.bam \
  --reference hg38/hg38.fa \
  --targets targets.bed
```

Use `-v` for high-level stage timings and GPU adapter selection. Use `-vv` when
you need more detailed per-input scanning diagnostics.

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

Common annotation options:

| Option | Purpose |
|---|---|
| `--database NAME=VCF` | Register an annotation database. |
| `--annotation NAME:SRC=DEST` | Copy an INFO field from a database into output records. |
| `--reference FASTA` | Normalize simple indel keys before matching. |
| `--index-type csi|tbi` | Select the output index type. |

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

Use `filter` after annotation or paired calling when the output VCF already
contains INFO or FORMAT fields needed for downstream triage.

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

Matching mode flags can be used with two-file and multi-set operations:

- default matching uses exact allele keys: `(CHROM, POS, REF, ALT)`
- `--site-only` matches by `CHROM/POS` and ignores alleles
- `--reference FASTA` left-normalizes non-symbolic allele keys before matching
- `--genotype-aware` keeps only ALT alleles supported by called non-reference
  `FORMAT/GT` values when samples are present

`--site-only` and `--reference` are mutually exclusive. Multi-ALT records are
split internally for matching; in genotype-aware mode, only GT-supported ALT
alleles contribute keys. Symbolic alleles, breakends, and spanning deletions are
not normalized by `--reference`. Filtered records are not treated specially; they
participate if their keys match. Missing genotypes such as `./.` or `.|.`,
missing alleles, and records without `FORMAT/GT` count as absence.

```bash
varlock intersect \
  --left caller_a.vcf.gz \
  --right caller_b.vcf.gz \
  --reference hg38.fa \
  --mode shared \
  --output normalized_shared.vcf.gz
```

```bash
varlock intersect \
  --left tumor.vcf.gz \
  --right panel.vcf.gz \
  --site-only \
  --mode left-only \
  --output novel_sites.vcf.gz
```

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

Multi-set mode compares named sets from repeated `--set` inputs or a tab-delimited
manifest. A set can use record presence for all records, or genotype support from
specific samples with `NAME=VCF:SAMPLE[,SAMPLE...]`.

```bash
varlock intersect \
  --set A=cohort1.vcf.gz:tumor_a,tumor_b \
  --set B=cohort2.vcf.gz:tumor_c \
  --set C=cohort3.vcf.gz:tumor_d \
  --mode all-shared \
  --output all_shared.vcf.gz
```

```bash
varlock intersect \
  --set A=file1.vcf.gz:a,b \
  --set B=file2.vcf.gz:c,d \
  --set C=file3.vcf.gz:e \
  --mode set-diff A-B \
  --output A_minus_B.vcf.gz
```

Manifest lines are `NAME<TAB>VCF[<TAB>SAMPLE[,SAMPLE...]]`:

```bash
varlock intersect \
  --set-manifest sets.tsv \
  --mode any-shared \
  --emit-set A \
  --output any_shared.vcf.gz
```

Multi-set output emits records from `--emit-set` or, by default, the first set
(`set-diff` defaults to the left side of `A-B`). Records include
`INFO/VARLOCK_SET_COUNT` and `INFO/VARLOCK_SETS`.

## Command Cheat Sheet

```bash
# Targeted calling
varlock call-targets -i sample.bam -r hg38.fa -T targets.bed -o calls.vcf.gz

# Force CPU in a GPU-enabled build
varlock call-targets --cpu -i sample.bam -r hg38.fa -T targets.bed

# List GPU adapters
varlock call-targets --gpu-list

# Paired tumor/normal calling
varlock call-targets -i tumor.bam -i normal.bam -r hg38.fa -T targets.bed \
  --pair tumor=Tumor,normal=Normal -o paired.vcf.gz

# Annotate from a VCF-like database
varlock annotate -i calls.vcf.gz --database db=db.vcf.gz \
  --annotation db:AF=db_AF -o annotated.vcf.gz

# Filter by INFO expression
varlock filter -i annotated.vcf.gz --expr "db_AF < 0.01 && DP >= 20" \
  -o rare.vcf.gz

# Intersect two VCFs
varlock intersect --left a.vcf.gz --right b.vcf.gz --mode shared \
  -o shared.vcf.gz
```

## Roadmap

Planned features are tracked with `bd` issues. This section lists known follow-up
areas rather than a complete release plan.

### Variant Intersection And Difference

Implemented:

- two-file exact-key `shared`, `left-only`, and `right-only`
- one-file multi-sample group `shared`, `left-only`, and `right-only`
- multi-file named-set `all-shared`, `any-shared`, and `set-diff A-B`
- matching modes: `--reference`, `--site-only`, and `--genotype-aware`

Follow-up:

- Add indexed or streamed multi-set evaluation for large cohorts.
- Decide whether projected VCF output or tabular summaries are needed for set
  operations.

### Annotation Follow-Ups

The annotation command supports exact and reference-normalized VCF-like database
matching, multi-ALT records, and indexed lookup for `.tbi` and `.csi` databases.
Future work should add TSV/BED-style variant or interval database adapters.

Design decisions to settle before broader database support:

- required database indexing and supported database formats
- exact allele matching vs site-only annotation
- optional normalization against a reference
- INFO/header naming and provenance lines
- no-hit representation
- multi-database conflict handling
