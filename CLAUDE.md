# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->


## Build & Test

```bash
cargo build                        # debug build
cargo build --release              # release build
cargo check                        # fast type-check without linking
cargo clippy                       # lint
cargo test                         # run all tests
cargo run -- call-targets --help   # run the CLI
```

No tests exist yet. The edition is 2024, so `let`-chains (`if let Some(x) = y && cond`) and other edition-2024 features are in active use throughout the codebase.

## Architecture Overview

`varlock` is a single-binary CLI for pileup-based SNV calling from BAM files. Currently only one subcommand exists: `call-targets`.

### Module responsibilities

**`main.rs`** — CLI layer only. Owns `CallTargetsArgs`, `ExecutionContext`, `IndexType`, `log_verbose`, and the Unix stderr-tee for `--log-file`. No bioinformatics logic lives here.

**`fasta_prep.rs`** — Self-contained reference preparation. `prepare_reference` is the sole public entry point: if the given FASTA already has a `.fai` (and `.gzi` for bgzipped inputs) it returns immediately; otherwise it sorts sequences by name, re-wraps at 60 bp, and writes a multi-threaded BGZF output alongside `.fai` and `.gzi` indices as `<name>.sorted.fa.gz` next to the original. Intermediate work is spooled to a temp dir that self-cleans on drop.

**`call_targets.rs`** — Currently a monolith that is intended to be split into focused sub-modules. Logical groupings within the file:

| Concern | Key items |
|---|---|
| Types / data model | `SiteKey`, `SiteCounts`, `Interval`, `TargetIndex`, `PreparedCallTargets` |
| Input resolution | `resolve_inputs`, `collect_bams_recursive` |
| Header / sample prep | `prepare_headers`, `collect_samples`, `rg_to_sm_from_header`, `read_rg_map` |
| Target loading | `load_targets` — parses BED, maps chroms to ref-seq IDs, merges overlapping intervals |
| BAM scanning / pileup | `process_input_bam`, `pileup_record`, `merge_counts`, `ScanParams`, `PileupSettings` |
| FASTA random access | `FastaIndex`, `open_fai`, `open_fasta_index` — line-cached seek into plain or BGZF-indexed FASTA |
| VCF output | `CallTargetsOutputState`, `begin/write_chunk/finish_call_targets_output`, `choose_alt` |

### Data flow for `call-targets`

```
CallTargetsArgs
  → fasta_prep::prepare_reference        (ensure .sorted.fa.gz + .fai + .gzi)
  → resolve_inputs                        (BAM paths from -i / --bamlist)
  → prepare_headers                       (validate shared ref seqs, build ref_name→id map)
  → load_targets                          (BED → TargetIndex with merged intervals)
  → collect_samples                       (RG→SM from BAM headers or --rg-map file)
  → for each BAM: process_input_bam       (pileup → BTreeMap<SiteKey, SiteCounts>)
  → merge_counts                          (accumulate across BAMs, cap at max_depth)
  → write_call_targets_output             (begin → chunk → finish: BGZF VCF + CSI/TBI index)
```

### Key invariants to preserve when refactoring

- `SiteKey.position` is **1-based** (VCF convention); BED intervals are **0-based half-open** — `in_targets` receives `ref_pos - 1` for the comparison.
- `SiteCounts.per_sample` is `[u32; 4]` indexed A=0, C=1, G=2, T=3 via `base_index`.
- `merge_sample_counts_with_cap` uses deterministic largest-remainder rounding when merging would exceed `max_depth`.
- VCF records must be written in `BTreeMap` iteration order (sorted by `(reference_sequence_id, position)`) — the CSI/TBI indexers require strictly increasing coordinates.
- The CSI/TBI indexer is updated per-record using BGZF virtual positions captured immediately before and after each `write_vcf_record` call.

## Conventions & Patterns

- All logging goes to stderr via `eprintln!` or `log_verbose`. Verbosity levels: `-v` = stage timings, `-vv` = per-file detail.
- `references/genemancer/` contains the upstream project this was ported from, including a GPU streaming path and additional subcommands. Consult it when adding features that were already implemented upstream.
