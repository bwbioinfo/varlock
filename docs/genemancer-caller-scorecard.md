# Genemancer Caller Capabilities — Varlock Adoption Scorecard

**Issue:** varlock-l4b.10  
**Status:** Provisional. Cargo is unavailable in this workspace; dynamic command and test execution cannot be performed. All findings are source-backed.  
**Frozen baseline reference:** `/workspace/genemancer/docs/scope-transfer-baseline.md` (2026-09-22)  
**Sources inspected (read-only):**
- `genemancer/src/main.rs`, `call_targets.rs`, `call_targets_gpu.rs`, `vcf.rs`, `Cargo.toml`
- `genemancer/src/call_targets_gpu/` submodules (aggregate_backends, aggregate_utils, progress, runtime, scan)
- `varlock/src/main.rs`, `call_targets/mod.rs`, `call_targets/gpu/mod.rs`, `call_targets/types.rs`, `call_targets/output.rs`, `call_targets/pileup.rs`, `call_targets/samples.rs`, `vcf.rs`, `Cargo.toml`
- `varlock/tests/call_targets_e2e.rs`, `tests/cli_smoke.rs`, `tests/support/`

---

## 1. Targeted Caller — CPU Path (`call-targets`)

### 1.1 CLI surface comparison

| Argument | Genemancer | Varlock | Notes |
|---|---|---|---|
| `--input` / `--bamlist` | `-i`/`--bamlist` (required unless other present) | `-i`/`--bamlist` (optional; no targets = all covered) | Varlock allows no inputs for no-target mode |
| `--reference` | Required `-r` | Optional `-r` (required at runtime) | Same semantics; Varlock enforces at runtime |
| `--targets` | Required `-T` | Optional `-T`; omit = whole genome | Varlock extends scope beyond targeted |
| `--output` | Required `-o` | Optional `-o`; default derived from first input | Varlock adds derived default |
| `--rg-map` | Optional RG→SM TSV | Optional RG→SM TSV; only valid with `--split-by sm` | Varlock adds explicit `--split-by` constraint |
| `--split-by` | Not present (always SM) | `sm` / `rg` / `file` | Varlock unique capability |
| `--index-type` | `csi` / `tbi`; default CSI | `csi` / `tbi`; default CSI | Identical |
| `--min-mapq` | 0–60, default 20 | 0–60, default 20 | Identical |
| `--min-baseq` | 0–60, default 20 | 0–60, default 20 | Identical |
| `--min-alt-count` | default 1 | default 1 | Identical |
| `--min-alt-fraction` | 0.0–1.0, default 0.0 | 0.0–1.0, default 0.0 | Identical |
| `--max-depth` | default 100,000 | default 100,000 | Identical |
| `--no-indels` / `--indels-only` | Not present (SNVs only) | Present | Varlock unique: indel support |
| `--pair` | Not present | `tumor=S,normal=S` with five sub-flags | Varlock unique: somatic paired calling |

**Summary:** The shared common subset (targeted SNV calling with BAM inputs, FASTA, BED, RG map, thresholds, indexed VCF.gz output) is structurally identical. Varlock has materially advanced the caller with: no-target whole-genome mode, three `--split-by` strategies, CIGAR-based indel calling, and paired somatic calling.

### 1.2 Internal architecture comparison

| Aspect | Genemancer | Varlock |
|---|---|---|
| CPU entry point | `call_targets::run(args, ctx)` delegates to `call_targets_gpu::run_cpu_streaming(args, ctx)` | `call_targets::run_cpu(args, ctx)` is self-contained in `call_targets/mod.rs` |
| Sample resolution | `collect_samples()` in `call_targets.rs`; single RG→SM strategy | `collect_sample_resolution()` with `InputSampleResolver` per input; supports SM/RG/file split strategies |
| Indel tracking | Absent from Genemancer CPU path | `IndelKey` + `IndelCounts` + CIGAR parsing in `pileup.rs` |
| Paired calling | Absent | `PairedCallingConfig` / `PairedSampleRoles` with full filter set |
| Output module | Monolithic in `call_targets.rs` (≈500 lines) | Extracted `call_targets/output.rs` (≈968 lines) with `OutputCall::Snv` / `OutputCall::Indel` merge-sort |
| Reference | `FastaIndex` with line cache inline in `call_targets.rs` | `call_targets/reference.rs` dedicated module |
| Depth policy open issue | Not gated (baseline notes P0 depth-policy work in Varlock) | `max_depth` cap; `merge_counts` with cap; no-target streaming path separate |

### 1.3 Functional completeness

Genemancer `call-targets`: targeted SNVs, CSI/TBI, RG→SM, MAPQ/BQ/alt thresholds, depth cap. Documented; implemented.

Varlock `call-targets`: superset — adds indels, no-target mode, split-by, paired calling, reference-normalized allele matching for annotation/intersect downstream. Documented CLI args confirmed in `main.rs`.

### 1.4 Correctness evidence

| Repo | Test mechanism | Coverage scope |
|---|---|---|
| Genemancer | No `call_targets`-specific integration tests found; only `cnloh_parity_baseline.rs` and `gtf_to_introns.rs` | CN-LOH fixtures and annotation unit tests only; targeted calling is untested at the integration level |
| Varlock | `tests/call_targets_e2e.rs` (15 test functions) + `tests/cli_smoke.rs` (≥5 call-targets smoke tests) = ≥20 `#[test]` attributes overall; e2e tests synthesize BAM fixtures, write FASTA/BED, invoke the binary, and assert VCF record fields, FORMAT/GT/DP/AD, CSI/TBI index structure | Basic SNV calling, RG/SM mapping, depth cap, index types, indel calling, `--split-by`, paired calling |

Varlock has demonstrably stronger fixture-based correctness evidence than Genemancer for this caller path. Genemancer has zero CPU `call-targets` integration tests.

**Verification limit:** `cargo test` was not executed. Test counts and assertions were confirmed by source inspection only.

### 1.5 Interchange semantics

Both repos use:
- 0-based, half-open BED coordinates as target intervals
- One-based VCF POS
- BGZF-compressed VCF output (`noodles_bgzf`)
- CSI or TBI index written immediately after the VCF
- RG→SM mapping via TSV

Varlock adds reference normalization via `vcf.rs::VariantNormalizer` (used in `annotate`/`intersect`). The Genemancer caller does not expose this, but does not contradict it.

### 1.6 Resource and operational behavior

| Aspect | Genemancer | Varlock |
|---|---|---|
| Thread model | Uses `--threads`/`ExecutionContext`; CPU streaming in `call_targets_gpu.rs` | Uses `--threads`/`ExecutionContext`; sequential BAM iteration in `run_cpu` |
| Memory bound | `max_depth` cap truncates per-site accumulation | Same `max_depth` cap; `merge_counts` applies cap |
| Streaming vs batch | Streaming scan in GPU module (`SCAN_BATCH_OBSERVATIONS = 100,000`) | Sequential scan over `BTreeMap<SiteKey, SiteCounts>` |
| Log/observability | `-v`/`-vv` stage timing via `log_verbose`; `--log-file` | Same flags and mechanism |

### 1.7 Provisional decision: `call-targets` CPU path

**RETIRE from Genemancer.** Varlock's CPU caller covers the full common targeted-SNV subset with identical thresholds and interchange semantics. Varlock adds indels, paired calling, and split-by strategies that are absent in Genemancer. Varlock has materially stronger test coverage (≥15 integration tests versus zero in Genemancer for this path). No duplicate algorithm maintenance is required.

**Gate:** The baseline notes that Varlock must resolve its P0 depth-policy and BAM/FASTA reference-validation work (varlock-l4b.2/.4) before Genemancer command deprecation. This evaluation does not change that gate. The depth-policy gap is internal to Varlock and does not block the retirement decision—it blocks the support-window start.

---

## 2. GPU/CUDA Surface (`call-targets-gpu`)

### 2.1 Genemancer GPU feature inventory

Source: `genemancer/src/call_targets_gpu.rs` (3,464 lines), `genemancer/src/main.rs` `CallTargetsGpuArgs`.

**Backend selection:** `--gpu-backend` with variants `auto`, `cuda`, `vulkan`, `metal`, `dx12`, `gl`, `browser-wgpu`; `--cuda-device INDEX`.

**Failure modes:** `--require-gpu` (fail if no GPU), default fallback to CPU.

**Diagnostics/tuning:**
- `--tuning-mode` (`auto`, `hybrid`)
- `--tuning-profile` (`balanced`, plus others)
- `--tuning-scale-percent` (1–500)
- `--wgpu-matrix-utilization-percent` (0–100)
- `--wgpu-upload-utilization-percent` (0–100)
- `--cuda-large-batch-mode` / `--no-cuda-large-batch-mode`
- `--defer-cuda-aggregation` / `--no-defer-cuda-aggregation`
- `--cuda-pileup-record-batch SIZE`
- `--cuda-pileup-base-budget SIZE`

**Experimental:**
- `--gpu-pileup` — experimental CUDA pileup extraction kernel (GPU-side CIGAR traversal). Enabled only with `--features cuda`.
- `--compare-cpu` — runs both GPU and CPU paths and reports mismatch counts; intended for validation only.

**CUDA kernel:** Inline PTX source in `call_targets_gpu.rs` — `aggregate_observations` kernel and `extract_observations` pileup kernel. Compiled at runtime via `cudarc` + `nvrtc`.

**wgpu kernel:** WGSL shader `AGGREGATE_SHADER` also inlined, for non-CUDA GPU backends.

**Streams and chunking:** Dedicated scan worker threads with `ScanEvent` channel; chunk planning based on matrix budget; deferred eager relief heuristics; dynamic tuning constants (≥25 named consts).

### 2.2 Varlock GPU feature inventory

Source: `varlock/src/call_targets/gpu/mod.rs` (620 lines), `runtime.rs`, `aggregate.rs`, `kernel.rs`, `scan.rs`.

**Feature gate:** `#[cfg(feature = "wgpu")]` throughout; CUDA is not present.

**Backend selection:** `--gpu-backend` with variants `all`, `vulkan`, `metal`, `dx12`, `gl` (no CUDA, no browser-wgpu).

**Failure modes:** `--require-gpu`, `--cpu` (force CPU path), `--gpu-list`.

**Adapter selection:** `--gpu-index` (by index), `--gpu-name` (by substring).

**Tuning:** `--matrix-budget-mib`, `--max-obs-upload`, `--obs-flush-threshold`.

**No-target GPU path:** Varlock adds `run_covered_gpu_path` for whole-genome streaming (absent in Genemancer).

**No CUDA, no experimental pileup:** Varlock's GPU module is wgpu-only; no CUDA feature, no inline PTX, no `--gpu-pileup`, no `--compare-cpu`.

### 2.3 Capability gap table

| Capability | Genemancer | Varlock | Gap direction |
|---|---|---|---|
| wgpu aggregate kernel (WGSL) | Yes | Yes | Parity |
| CUDA backend (`cudarc`, `nvrtc`) | Yes (`--features cuda`) | No | Genemancer ahead |
| Experimental GPU pileup (`--gpu-pileup`) | Yes (CUDA only) | No | Genemancer ahead (experimental) |
| CPU–GPU comparison mode (`--compare-cpu`) | Yes | No | Genemancer ahead |
| No-target (whole-genome) GPU path | No | Yes (`run_covered_gpu_path`) | Varlock ahead |
| Backend: auto/cuda/vulkan/metal/dx12/gl | Yes (`--gpu-backend`) | Subset: all/vulkan/metal/dx12/gl | Genemancer has CUDA backend variant |
| Per-adapter index/name selection | No | Yes (`--gpu-index`, `--gpu-name`) | Varlock ahead |
| Tuning knobs (matrix/upload budgets) | Verbose set (8+ flags) | Simplified set (3 flags) | Genemancer broader; Varlock cleaner |
| GPU correctness test | None found | `call_targets_e2e.rs`: `#[cfg(feature = "wgpu")]` block asserts GPU alleles equal CPU alleles | Varlock ahead |
| Tier classification (Datacenter/HighEnd/Discrete/Integrated) | Not found in Genemancer source | `GpuTier` in `runtime.rs` with auto-tuning per tier | Varlock ahead |

### 2.4 CUDA surface assessment

The Genemancer `--features cuda` surface is experimental and non-portable:

- The `--gpu-pileup` flag is described in `scope-transfer-baseline.md` as experimental.
- CUDA compilation is runtime (nvrtc), creating a deployment dependency on `libcuda` and `libnvrtc`.
- The `--compare-cpu` flag is a developer/validation tool, not a user feature.
- No CUDA integration tests were found in Genemancer.
- The `cudarc` dependency version (`0.19.3`) pinned to `cuda-12080` with fallback-dynamic-loading is a non-trivial maintenance surface.

**The CUDA-specific capabilities do not constitute a material, verified, safely portable advantage** over Varlock's wgpu path under the rubric criteria (correctness evidence, interchange semantics, maintenance cost). Genemancer's wgpu aggregate kernel is functionally equivalent to Varlock's; the only Genemancer-unique stable GPU contribution is the broader `--gpu-backend` enum (CUDA + browser-wgpu variants).

### 2.5 Provisional decision: `call-targets-gpu`

**RETIRE from Genemancer.** Varlock's `call-targets-gpu` (wgpu feature) is the safe destination. The common wgpu aggregate kernel, failure-mode semantics (`--require-gpu`/fallback), and SNV output contract are already present in Varlock.

**CUDA surface: document as unsupported for migration.** The CUDA backend, `--gpu-pileup`, and `--compare-cpu` are Genemancer-only. They must be:
1. Documented explicitly in the migration guide as not carried to Varlock.
2. Not replicated in Varlock unless a separate tracked evaluation confirms correctness and maintainability.
3. Considered for potential future adoption only after Varlock's wgpu path is validated at P0 correctness level.

**Varlock gaps to fill before retirement:** `--gpu-backend` in Varlock lacks `cuda` and `browser-wgpu` variants. The missing variants do not block migration for CPU/GPU equivalent SNV output, but the migration guide must document them as unsupported. The no-target GPU path in Varlock is a Varlock advantage, not a migration blocker.

---

## 3. VCF Operations (`vcf diff`)

### 3.1 Genemancer `vcf diff` state

Source: `genemancer/src/vcf.rs`.

The `vcf diff` command validates input semantics (multisample vs. multi-file mode, set definitions, member uniqueness, basename disambiguation) and then explicitly fails:

```
bail!("vcf diff scaffold: input semantics are validated, but record loading and set-difference computation are not implemented yet")
```

There are no VCF records read, no diff computation, and no output produced. The input-validation logic has unit tests (6 tests in `vcf.rs` `mod tests`).

### 3.2 Varlock VCF operations

Varlock provides `intersect`, `filter`, and `annotate`:

- **`intersect`**: Two-file mode (`--left`/`--right`), one-file multi-sample mode (`--left-samples`/`--right-samples`), multi-set mode (`--set`/`--set-manifest`/`--emit-set`), modes `shared`/`left-only`/`right-only`/`all-shared`/`any-shared`/`set-diff`, `--genotype-aware`, `--site-only`, reference normalization. Smoke tests confirm set-diff produces correct output.
- **`filter`**: INFO predicates (`--require-info`, `--exclude-info`, `--max-info`, `--expr`), sample/group genotype filters. Smoke tested.
- **`annotate`**: Database annotation with INFO field mapping. Smoke tested with BGZF + CSI output verification.

### 3.3 Provisional decision: `vcf diff`

**RETIRE from Genemancer.** The command is a scaffold with no record-processing result. There is no output parity obligation per the baseline. Varlock's `intersect --mode set-diff` provides the functional destination. The input-validation error contract may be preserved where useful (i.e., clear error messages for invalid set definitions), but no code porting is required or warranted.

---

## 4. Test and Validation Maturity Summary

| Component | Genemancer tests | Varlock tests | Verdict |
|---|---|---|---|
| `call-targets` CPU | 0 integration tests | ≥15 integration tests (e2e): SNV, RG/SM, depth, index, indels, split-by, paired | Varlock materially stronger |
| `call-targets-gpu` | No GPU integration tests found | `#[cfg(feature = "wgpu")]` test block asserts GPU/CPU output parity | Varlock ahead |
| VCF operations | 6 unit tests (input validation only; no record tests) | CLI smoke tests for annotate, filter, intersect (with record-level assertions) | Varlock stronger; Genemancer has no record-level VCF tests |
| Total `#[test]` in scope | 3 (cnloh parity baseline + 71 in whole repo) | ~20 in test files + ~10 inline | Genemancer test volume is dominated by CN-LOH, not caller paths |

**Verification limit:** Test counts from source line counts and `grep` on `#[test]` attributes. No `cargo test` run was performed.

---

## 5. Resource, Reference, and Sample Semantics

### Reference handling

Both repos use `fasta_prep::prepare_reference()` with the same signature and semantics: accepts plain FASTA or BGZF, creates `.fai` if absent, returns a prepared path. The modules are structurally identical (both `src/fasta_prep.rs` present in both repos).

### Sample identity

- Genemancer: SM tag from BAM header RG records, optional override via `--rg-map` TSV (RG\tSM). One strategy.
- Varlock: Same SM-tag default, plus `--split-by rg` (one column per RG) and `--split-by file` (one column per input BAM), with `InputSampleResolver` abstraction per input. Superset.

### Coordinate convention

Both: 0-based, half-open BED targets. Both: 1-based VCF POS. Both use `noodles_core::Position` for VCF record writing.

### Read filtering

Both: `--min-mapq` (MAPQ filter), secondary/supplementary/unmapped flags skipped. Confirmed by pileup source inspection.

---

## 6. Maintenance Cost

| Factor | Assessment |
|---|---|
| Duplicate algorithm risk | Genemancer `call_targets_gpu.rs` is 3,464 lines; `call_targets.rs` is monolithic. Varlock splits functionality into 8+ focused modules. Maintaining both would require synchronizing threshold logic, RG resolution, output format, and index generation across both repos. |
| CUDA dependency | `cudarc` 0.19.3 locked to `cuda-12080` is a significant maintenance surface. Varlock does not carry this. |
| Unique Genemancer capabilities | `--gpu-pileup` (experimental, no tests), `--compare-cpu` (developer tool, no tests), CUDA backend. None have production correctness evidence. |
| Varlock unique capabilities | Indels, no-target mode, split-by, paired calling, no-target GPU path, per-adapter selection, tier auto-tuning. All are production-oriented with test coverage. |

**Adopting from Genemancer into Varlock is not recommended for any component.** Varlock's implementation already supersedes the common caller subset. The Genemancer-only capabilities are experimental and untested.

---

## 7. Per-Component Decision Table

| Genemancer surface | Decision | Rationale |
|---|---|---|
| `call-targets` (CPU) | **RETIRE** | Varlock superset; zero Genemancer integration tests; identical interchange semantics. Gate: varlock-l4b.2/.4 must close first. |
| `call-targets-gpu` (wgpu) | **RETIRE** | Varlock wgpu path covers the common SNV subset with GPU correctness test. |
| `call-targets-gpu` — CUDA backend | **UNSUPPORTED** | Experimental, no correctness evidence, non-portable runtime dependency. Document in migration guide. |
| `call-targets-gpu` — `--gpu-pileup` | **UNSUPPORTED** | Experimental CUDA feature, no tests. Do not carry to Varlock. |
| `call-targets-gpu` — `--compare-cpu` | **UNSUPPORTED** | Developer/validation tool, not a user feature. Do not carry; equivalent manual comparison possible. |
| `vcf diff` | **RETIRE** | Scaffold only; no record processing; Varlock `intersect --mode set-diff` is the destination. |

---

## 8. Open Blockers and Follow-up

1. **varlock-l4b.2/.4 (P0 depth-policy / reference-validation):** The baseline explicitly requires these to close before Genemancer `call-targets` is deprecated. This scorecard does not change that gate.
2. **Migration guide:** A command-option mapping document is required before the support window begins. It must explicitly list `--gpu-pileup`, `--compare-cpu`, and `--gpu-backend cuda` as unsupported in Varlock.
3. **Dynamic CLI verification:** Because Cargo is unavailable, `--help` output, version strings, and fixture checksums could not be captured. These must be run in a Rust-enabled environment per the reproducible baseline commands in `scope-transfer-baseline.md` before the migration is executed.
4. **`call-targets-gpu` no-target parity:** Varlock has `run_covered_gpu_path`; Genemancer does not. This is a Varlock advantage, not a blocker, but should be confirmed when the whole-genome mode is exercised in Varlock integration tests.

---

*Scorecard created 2026-09-22. Source-backed only; no binaries were executed.*
