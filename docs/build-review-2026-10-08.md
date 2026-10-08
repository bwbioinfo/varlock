# Varlock build review, October 8, 2026

## Scope and verdict

Reviewed checkpoint `61f42b3b720855a72719f11e2ddfe89c8ad0a00e` on the isolated Jbox checkout. The checkout was clean at review start. NanoCov work is shelved.

**The reported default build failure is reproducible and is a source/API error, not merely missing dependencies.** The optional GPU build has an additional compiler-version blocker. No production Rust source, dependency versions, lockfile or workflow was changed during this review. The reader repair described below was tested only in a scratch copy.

## Findings, in priority order

### 1. P1: Invalid reader trait object prevents every default binary build

Location: `src/filter/native.rs:52-55`, introduced by checkpoint `61f42b3`.

`cargo build --release --locked` fails with E0107 at the `Box<dyn vcf::variant::io::Read>` declaration. Locked `noodles-vcf 0.85.0` defines `Read<R>`, not a nongeneric `Read`. Both VCF and BCF readers implement that trait for their own buffered transport type. The native filter module is compiled unconditionally via `src/filter.rs:1`, so this prevents the entire CLI from compiling, even when no BCF operation is requested.

Merely adding a generic argument without reconciling the two transports is insufficient. A minimal repair can erase both transports to `Box<dyn BufRead>`, use `Read<Box<dyn BufRead>>`, and construct the BCF reader with `Reader::from` around an explicitly BGZF-decoding boxed transport. Alternatively, an enum can dispatch the two concrete readers. Do not replace BGZF decoding with a raw file reader.

**Validation of the proposed repair:** An isolated copy with only that transport/type repair compiled in release mode and passed all **125 default tests**, including all seven native BCF end-to-end tests. This does not mean the original checkout is repaired.

### 2. P1 for GPU builds here: WGPU requires a newer compiler than the guest provides

Location: `Cargo.toml:33`, `Cargo.lock:1502-1504`; toolchain requirements absent from `Cargo.toml` and nonspecific in `README.md:75-81`.

The optional dependency is locked to `wgpu 28.0.0`. Its downloaded manifest declares `rust-version = "1.92"`. The guest provides `rustc 1.89.0` and has no `rustup` command. Actual `cargo check --locked --all-targets --features wgpu` fails at Cargo's compiler compatibility check before Varlock's GPU source is type-checked.

Document and provision Rust 1.92 or newer for this feature, or deliberately select a compatible WGPU/API version and verify it. Do not blindly downgrade the lockfile or use `--ignore-rust-version` as a claimed fix. A newer compiler alone cannot fix finding 1. GPU source compilation and GPU runtime behavior remain unverified in this guest.

### 3. P2: The manual all-feature CI job cannot run

Location: `.github/workflows/ci.yml:3-6` and `:43-48`.

The workflow declares only push and pull-request triggers. Its all-feature job requires `github.event_name == 'workflow_dispatch'`, an event the workflow never enables. Consequently, neither normal CI nor the intended manual trigger can execute that job. The comment about special GPU runners does not make the trigger reachable.

Add the intended manual trigger or change the scheduling rule. Separately, a compile-only GPU-feature gate can run without exercising a GPU adapter. Keep runtime/device checks distinct from compile checks.

### 4. P2: Strict quality gates remain broken after the reader repair

Locations: `src/call_targets/pileup.rs:18`, `tests/filter_bcf_e2e.rs:92`, and the formatting gate `.github/workflows/ci.yml:21`.

In the reader-repaired isolated copy, strict all-target Clippy still fails on an unread `PileupSettings.max_depth` field and an unnecessary mutable `query` binding in the BCF end-to-end test. CI also sets `RUSTFLAGS="-D warnings"`, so these warnings cannot be treated as harmless under the configured policy. Removing the unused field should not remove the actual final merged-count depth cap.

The untouched checkout also fails `cargo fmt --all -- --check` in eight distinct source/test files. Restore the configured gates rather than lowering their strictness. These are independent CI blockers, not the root cause of ordinary non-strict release compilation.

## Actual evidence

| Check | Result | Scope |
| --- | --- | --- |
| `cargo build --release --locked --offline` | Cache missing `noodles-bcf` | Original checkout, before restore |
| `cargo build --release --locked` | E0107 at native reader declaration | Original checkout, locked dependencies restored |
| `cargo check --locked --all-targets --features wgpu` | WGPU requires Rust 1.92, installed 1.89 | Original checkout |
| `cargo fmt --all -- --check` | Failed in eight distinct files | Original checkout |
| Release build with minimal reader repair | Passed | Isolated scratch copy only |
| Default tests with minimal reader repair | 125 passed, zero failed/ignored | Isolated scratch copy only |
| Strict all-target Clippy after reader repair | Failed on dead field and unused mut | Isolated scratch copy only |

Independent read-only review corroborated the source API mismatch, WGPU compiler requirement and unreachable all-feature CI trigger. Logs are retained under `/home/jbox/.jcode/scratch/` with names `varlock-release-baseline-20261008.log`, `varlock-release-build-20261008.log`, `varlock-gpu-check-20261008.log`, `varlock-fmt-review-20261008.log`, and `varlock-isolated-reader-probe-20261008.log`.

The isolated probe is `/home/jbox/.jcode/scratch/varlock-review-probe-20261008`. Its binaries and tests are not evidence that the unmodified production checkout builds.

## Reproduce the principal blockers

The following commands are standalone and target the original checkout, not the repaired probe:

```bash
cd /workspace/varlock && CARGO_BUILD_JOBS=2 cargo build --release --locked
```

```bash
cd /workspace/varlock && CARGO_BUILD_JOBS=2 cargo check --locked --all-targets --features wgpu
```

```bash
cd /workspace/varlock && cargo fmt --all -- --check
```

## Recommended next step

Apply the small reader transport/type repair to the real checkout, restore strict lint/format gates, and rerun default tests. Then explicitly provision the GPU-supported compiler and verify both GPU-feature compilation and relevant runtime paths. Track integration and real-data correctness separately from this build diagnosis.

## Toolchain follow-up, October 8, 2026, 22:10 UTC

The user requested Rust 1.99.0 in the shared workspace Dockerfile and delegated
GPU runtime tests to their environment outside the container. The Dockerfile
now pins `rust:1.99.0-bookworm`. Rust/Cargo 1.99.0 with matching rustfmt and Clippy
are also installed in executable workspace storage for this existing guest.
Select them with `source /workspace/project/.jbox/rust-env.sh`.

On Rust 1.99.0, the original checkout's all-target GPU-feature check reaches
the same reader E0107 instead of failing the compiler-version requirement.
The isolated reader-repaired copy passes `cargo check --locked --offline
--all-targets --features wgpu`. No GPU runtime test was executed and no Varlock
production source was repaired. The log is
`/home/jbox/.jcode/scratch/varlock-rust199-compile-20261008.log`.
The actual image rebuild still requires host Docker/Jbox, which is unavailable
inside this guest.
