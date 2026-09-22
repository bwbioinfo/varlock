# Varlock VCF and BCF processing contract

**Contract version:** 1.0 (2026-09-22)
**Status:** VCF contract retained. A native `filter` BCF implementation is present
in source but has not been compiled or runtime-validated in this environment.
BCF is not yet promoted to the supported contract.

This document is the user-facing capability contract for Varlock's `call-targets`,
`annotate`, `filter`, and `intersect` commands. It is based on the current CLI,
Rust implementation, tests, and the Genemancer adoption scorecard. The matrix
retains the VCF contract. The separately labeled implementation note records
source changes without claiming unexecuted tests passed.

## Executive summary

- Varlock currently reads and writes **text VCF**. It accepts plain text VCF and
  gzip-family input selected by the filename extension `.gz`, `.bgz`, or `.bgzf`.
- Varlock's generated variant output is **BGZF-compressed VCF**, accompanied by
  either a CSI index (default) or a TBI index (`--index-type tbi`). Output paths
  are documented as `VCF_GZ`; Varlock does not infer or emit BCF from an output
  suffix.
- **BCF is not yet a validated all-command capability.** `filter` now has a
  native noodles BCF source path selected by explicit format flags, described
  below. `call-targets`, `annotate`, and `intersect` remain VCF-only. A binary
  BCF file is never valid input to the legacy line-oriented VCF reader.
- No command preserves arbitrary input header bytes. `filter` and `intersect`
  copy header lines as text while selecting records; `annotate` copies header
  lines and adds annotation declarations as needed; `call-targets` creates a new
  Varlock header. Records are emitted as VCF text, with the documented
  transformations below.

## Supported contract matrix (not expanded by unexecuted BCF tests)

| Command | Plain VCF input | gzip-family VCF input | BCF input | Output | Index | Streaming/resource behavior |
|---|---:|---:|---:|---|---|---|
| `call-targets` | N/A | N/A | N/A | BGZF VCF | CSI default, or TBI | Sequential BAM/pileup processing; caller applies `--max-depth`; output is streamed through BGZF and indexed while written |
| `annotate` | Yes | Yes, when extension is `.gz`, `.bgz`, or `.bgzf` | No | BGZF VCF | CSI default, or TBI | Reads the input sequentially; indexed BGZF VCF/CSI/TBI databases can be queried without a reference normalizer; normalized lookup may load database annotations in memory |
| `filter` | Yes | Yes, when extension is `.gz`, `.bgz`, or `.bgzf` | No | BGZF VCF | CSI default, or TBI | Sequential record processing; records are not all loaded by the filter path |
| `intersect` | Yes | Yes, when extension is `.gz`, `.bgz`, or `.bgzf` | No | BGZF VCF | CSI default, or TBI | Two-file and set operations load variant keys/sets; selected output records are streamed; memory is proportional to the key/set data, not the emitted output |

The `.gz` reader uses a multi-member gzip decoder. The output writer is
`noodles-bgzf`, so output is BGZF rather than an arbitrary gzip stream.
Filename suffixes do not make a binary BCF file valid VCF input.

## Native filter implementation note (source-reviewed, runtime unvalidated)

This is an implementation inventory, **not a claim of passing runtime validation
or a promotion of contract version 1.0**. Cargo, rustc, and rustfmt were unavailable
when this slice was authored. `varlock-mku.1.1` tracks compilation, formatting,
testing, interoperability, and any resulting fixes.

- `filter --input-format vcf|bcf --output-format vcf|bcf` selects explicit codecs.
  Both default to `vcf`, retaining the existing text VCF path when neither
  endpoint is BCF. A `.bcf` filename alone never selects BCF.
- The new path uses `noodles-bcf` 0.83 and `noodles-vcf` 0.85 variant interfaces.
  BCF input uses a BGZF reader regardless of suffix. Raw uncompressed BCF and
  ordinary gzip BCF are not offered as input modes. VCF input retains the
  existing plain/gzip extension policy. Both output encodings are BGZF.
- Processing is sequential and requires no input index. One record is rendered
  in memory as VCF for existing INFO/sample predicates. Selected BCF output is
  written from the native record, not reparsed predicate text. There is no shell
  or temporary-file conversion. VCF output is the noodles rendering, not original
  text bytes.
- Headers are parsed and reserialized. Source contig, INFO, FORMAT, and sample
  definitions supply the native dictionaries. Missing declarations are not
  invented for BCF output. Numeric widths and missing values obey the codecs.
  Exact numeric spelling, header ordering, and binary byte identity are not
  promised.
- BCF output rejects `--index-type tbi` before output creation. It writes
  `<output>.csi` with native reference IDs, BGZF chunks, and variant start/end
  spans. This CSI has no Tabix text header (min_shift 14, depth 6). Selected BCF
  records must be sorted by output dictionary ID and position, with a usable
  position and end. VCF output retains the existing CSI/TBI writer.
- Native-path output and both possible index paths must not already exist,
  preventing input alias truncation and stale sidecars. Header, option, and
  sample validation precede data-file creation. Later record/output/index errors
  may leave partial files. Publication is not transactional. Do not consume
  failed command output.
- `tests/filter_bcf_e2e.rs` adds conversion, BCF magic, multiallelic and missing/
  phased genotype, sample/INFO filtering, empty-selection, native CSI span query,
  and rejection assertions for TBI, wrong codec, stale output, missing contigs,
  unsorted records, truncated records, and other-command format flags. These
  assertions have **not been executed** here.

`varlock-mku.1.2` tracks native BCF output for `call-targets` (including split/paired
outputs), primary/database BCF input and output for `annotate` with dictionary
remapping and database lookup, and BCF input/output across all `intersect` modes.
These commands do not accept the filter-only format flags and their VCF paths
are unchanged. Parent `varlock-mku.1` remains incomplete until runtime validation
and all-command acceptance criteria are met.

## Compression and indexing (existing VCF contract)

### Output

All four workflows that produce variants write BGZF-compressed VCF and then
write the requested sidecar index. The default is CSI. `--index-type tbi`
selects a Tabix VCF index. The sidecar names are formed by appending `.csi` or
`.tbi` to the complete output filename, for example:

```text
calls.vcf.gz       # BGZF VCF
calls.vcf.gz.csi   # default
# or calls.vcf.gz.tbi with --index-type tbi
```

Output records must be coordinate-addressable VCF records with a valid CHROM
and one-based POS for index construction. A malformed record or an index write
failure is an error; Varlock does not silently produce an unindexed partial
contract output.

### Input indexes

`annotate` can use a `.tbi` or `.csi` sidecar for a BGZF-compressed annotation
VCF when no reference normalizer is requested. The implementation looks for
`<input>.tbi` and `<input>.csi`. A missing database index falls back to an
in-memory read path where supported. `filter` and `intersect` do not require an
input index because they scan or materialize their own key sets. `annotate` does
not consume BCF database indexes. The new filter source path scans BCF without
an input index, as described above.

## Headers, source records, samples, and genotypes

### Header behavior

- `call-targets` writes a new VCFv4.3 header containing `##source=varlock
  call-targets`, the requested reference, contig declarations, Varlock INFO and
  FORMAT declarations, and a generated `#CHROM` line. It does not copy a source
  VCF header because its source is BAM plus FASTA.
- `filter` copies input header lines and emits selected record lines. It does
  not rewrite the header to describe filter expressions.
- `intersect` emits the header from the selected source set and emits selected
  source records. It does not merge arbitrary metadata from all input headers.
- `annotate` copies the input header and adds declarations for annotation fields
  written by the command. It does not preserve byte-for-byte header ordering as
  a compatibility guarantee.

Header lines and records are treated as UTF-8 text lines. Header metadata that
is not explicitly changed remains semantically present where the command copies
the input header, but exact whitespace, ordering, and byte identity are not a
stable API.

### Record and sample behavior

`filter` and `intersect` preserve selected source record text, except for the
normal output re-encoding through BGZF. `annotate` preserves source record
fields and updates INFO values according to the requested mappings. Existing
sample columns and FORMAT/sample values are retained by these text paths.

The filter and intersect predicates interpret samples by the `#CHROM` sample
order and FORMAT keys. Unknown sample names, malformed option expressions,
records with fewer than eight VCF fields, invalid positions, and inconsistent
sample references are errors rather than guessed behavior. Genotype-aware
intersection requires alternate support in the available genotype fields; it
is not a general genotype normalization engine.

`call-targets` creates samples from its BAM/read-group/sample-resolution policy.
Its generated records use the caller's documented `GT:DP:AD` semantics. SNVs
are selected from non-reference support; observed indels are emitted as
biallelic records. Symbolic alleles and structural variants are not emitted by
this caller path.

## Coordinates and reference normalization

- VCF `POS` is one-based. BED target intervals are zero-based, half-open.
- Without `--reference`, annotation and intersection use literal variant keys:
  CHROM, POS, REF, and each ALT allele. Matching does not rename contigs or
  normalize alleles.
- With `--reference`, `annotate` and `intersect` use the FASTA-backed
  `VariantNormalizer`. It uppercases ordinary alleles, trims common suffixes
  and prefixes, and left-normalizes indels using the reference sequence. Symbolic
  or special alleles are left unchanged. The normalized key affects matching,
  not the source record text written to output.
- The reference must be available through Varlock's indexed FASTA access path,
  and contig names must match the records being normalized. Varlock does not
  promise automatic `chr` prefix conversion, liftover, or reference assembly
  detection.
- `call-targets` anchors simple indels using the preceding aligned reference
  base, but its emitted calls are not left-normalized. This distinction is
  intentional and must not be treated as an equivalence guarantee with
  normalized annotation/intersection keys.

## Errors and unsupported records

Varlock fails with a contextual error for unreadable files, malformed VCF
lines, missing required headers, invalid coordinates, invalid options, missing
referenced samples, incompatible indexes, failed BGZF output, and failed index
construction. It does not promise recovery from malformed VCF or conversion of
BCF by inspecting its filename.

The current text implementation is not a general VCF validator. Unsupported
or only partially modeled constructs in the VCF-only commands include BCF, symbolic/structural
calling output, arbitrary binary encodings, and transformations requiring full
schema-aware rewriting of every INFO/FORMAT cardinality. Callers should
validate complex records with a dedicated VCF/BCF validator before relying on
semantic preservation.

## Provenance and compatibility policy

Generated caller output identifies itself with `##source=varlock call-targets`
and records the reference path in `##reference`. Transforming commands retain
source headers where their current implementation copies them, but do not add a
stable provenance chain for every input or command option. Therefore consumers
must treat the command line, Varlock version, and this contract version as the
complete provenance record for a run.

This is contract version 1.0. Patch releases may fix parsing, indexing, or error
handling without changing declared format support. A minor contract revision may
add a supported VCF feature or an explicitly documented normalization while
retaining existing valid workflows. A major revision is required to add BCF,
change default compression/index behavior, change coordinate semantics, or
change source-record/header preservation in a way that can alter downstream
results. New behavior must be covered by conversion/rejection tests before it
is advertised as supported.

## Migration from Genemancer

Retire Genemancer's VCF-diff/VCF-processing path only for workflows covered by
this contract. Invoke Varlock's `annotate`, `filter`, or `intersect` on plain or
BGZF VCF and consume BGZF VCF plus CSI/TBI output. Preserve Genemancer's
one-based VCF and zero-based half-open BED conventions. If a workflow currently
starts or ends in BCF, add an explicit external conversion boundary, for
example with a format-aware tool, and validate the converted VCF before and
after Varlock. Do not claim a drop-in BCF migration until the follow-up BCF
support work is complete and tested.

For Genemancer caller migration, Varlock's targeted CPU path provides the
common BAM/FASTA/BED caller behavior and indexed BGZF VCF output documented in
the adoption scorecard. Varlock additionally supports indels, whole-covered
position mode, split-by policies, and paired calling. Caller output is newly
constructed, so Genemancer source VCF headers are not applicable to that path.

## Known implementation gaps

1. **BCF coverage is incomplete and unvalidated.** A native filter reader/writer
   and focused tests are present in source. Runtime validation is tracked by
   `varlock-mku.1.1`; the three remaining command paths by `varlock-mku.1.2`.
   Source presence and unexecuted tests do not establish supported behavior.
2. **No end-to-end contract fixture covers all VCF input/output combinations.**
   Existing source and unit/e2e tests cover selected VCF, BGZF, CSI, and TBI
   cases, but do not yet provide a single matrix proving every command's
   accepted/rejected format behavior.
3. **Provenance is limited.** Transforming commands do not record a complete
   command-line/provenance chain in output metadata.

Follow-up Beads are filed for the first two gaps. Provenance is documented as a
policy limitation and should not be silently expanded in this contract.

## Source basis

- `src/main.rs`: command arguments, plain/gzip VCF wording, output/index options.
- `src/vcf.rs`: extension-based text reader, variant keys, normalizer, index
  construction, and malformed-record errors.
- `src/annotate.rs`, `src/filter.rs`, `src/intersect.rs`: command-specific
  reading, matching, filtering, annotation, and output paths.
- `src/call_targets/output.rs`: BGZF writer, generated header, caller records,
  and CSI/TBI sidecars.
- `src/filter/native.rs`, `tests/filter_bcf_e2e.rs`: native filter implementation
  and unexecuted conversion/rejection assertions.
- `Cargo.toml`, `Cargo.lock`: compatible noodles BCF/VCF production dependencies.
  The added BCF lock entry was checked against published metadata, not generated
  or resolved by Cargo in this environment.
- `README.md` and `docs/genemancer-caller-scorecard.md`: public behavior,
  coordinate conventions, and Genemancer migration baseline.
