#!/usr/bin/env bash
# local-ci.sh — Replay the full CI suite locally, including all-features tests.
#
# Usage:
#   ./local-ci.sh [--sha <SHA>] [--allow-dirty]
#   --sha must equal checked-out HEAD (assertion only). Dirty tree is refused
#   unless --allow-dirty, which forces overall=fail (never sha-bound evidence).
#
# Requirements: Rust stable toolchain with rustfmt + clippy components.
#   rustup component add rustfmt clippy
#
# Exit code: 0 = all gates passed, non-zero = at least one gate failed.
#
# Evidence record is written to local-ci-evidence.json at the repo root.
# Clean runs are SHA-bound; dirty diagnostic runs are explicitly marked unbound.

set -euo pipefail

# Scope Git checks and evidence to this script's repository, even when invoked
# from a subdirectory or another working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

REQ_SHA=""
ALLOW_DIRTY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sha) [[ $# -ge 2 && -n "$2" ]] || { echo "Usage: $0 [--sha <SHA>] [--allow-dirty]" >&2; exit 2; }
           REQ_SHA="$2"; shift 2 ;;
    --allow-dirty) ALLOW_DIRTY=1; shift ;;
    *) echo "Usage: $0 [--sha <SHA>] [--allow-dirty]" >&2; exit 2 ;;
  esac
done

# Evidence is bound to the checkout actually tested. --sha is only an
# assertion: it must resolve to the checked-out HEAD, never to another commit.
SHA="$(git rev-parse HEAD 2>/dev/null)" || { echo "error: not a git checkout" >&2; exit 2; }
if [[ -n "$REQ_SHA" ]]; then
  REQ_FULL="$(git rev-parse --verify --quiet "${REQ_SHA}^{commit}" 2>/dev/null || true)"
  if [[ "$REQ_FULL" != "$SHA" ]]; then
    echo "error: --sha $REQ_SHA does not match checked-out HEAD $SHA; check it out first" >&2
    exit 2
  fi
fi
# Dirty = tracked changes or untracked files across the entire repository.
dirty_state() { git status --porcelain --untracked-files=all -- . ':(exclude)local-ci-evidence.json' | wc -l | tr -d ' '; }
DIRTY_START="$(dirty_state)"
if [[ "$DIRTY_START" != "0" && $ALLOW_DIRTY -ne 1 ]]; then
  echo "error: working tree is dirty ($DIRTY_START paths); commit/stash or pass --allow-dirty (evidence will be marked dirty)" >&2
  exit 2
fi


EVIDENCE_FILE="local-ci-evidence.json"
START_TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
RESULTS=()

run_gate() {
  local name="$1"; shift
  local t0; t0=$(date +%s%3N)
  local exit_code=0
  local log
  log=$(mktemp)
  "$@" >"$log" 2>&1 || exit_code=$?
  local t1; t1=$(date +%s%3N)
  local dur=$(( t1 - t0 ))
  local status="pass"
  [[ $exit_code -ne 0 ]] && status="fail"
  printf '[%s] %-30s %s  (%dms)\n' "$status" "$name" "(exit $exit_code)" "$dur"
  # Store as a JSON fragment (no jq required)
  RESULTS+=("{\"gate\":\"$name\",\"status\":\"$status\",\"exit_code\":$exit_code,\"duration_ms\":$dur}")
  if [[ $exit_code -ne 0 ]]; then
    echo "--- STDOUT/STDERR ---"
    cat "$log"
    echo "--- END ---"
  fi
  rm -f "$log"
  return $exit_code
}

echo "=== varlock local CI ==="
echo "SHA: $SHA (dirty paths: $DIRTY_START)"
echo "Started: $START_TS"
echo ""

OVERALL=0
run_gate "rustfmt"          cargo fmt --all -- --check          || OVERALL=1
run_gate "clippy"           cargo clippy --all-targets -- -D warnings || OVERALL=1
run_gate "test-default"     cargo test                          || OVERALL=1
run_gate "test-all-features" cargo test --all-features          || OVERALL=1

END_TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
# Evidence is sha-bound only if HEAD and the tree stayed clean throughout.
HEAD_END="$(git rev-parse HEAD)"
DIRTY_END="$(dirty_state)"
VERIFIED=true
if [[ "$HEAD_END" != "$SHA" || "$DIRTY_START" != "0" || "$DIRTY_END" != "0" ]]; then
  VERIFIED=false
  OVERALL=1
  echo "warning: checkout changed or was dirty; evidence is NOT sha-bound, forcing FAIL" >&2
fi

# Write evidence JSON
{
  printf '{\n'
  printf '  "schema": 1,\n'
  printf '  "sha": "%s",\n' "$SHA"
  printf '  "sha_bound": %s,\n' "$VERIFIED"
  printf '  "dirty_paths_start": %s,\n' "$DIRTY_START"
  printf '  "dirty_paths_end": %s,\n' "$DIRTY_END"
  printf '  "environment": "local",\n'
  printf '  "hosted": false,\n'
  printf '  "note": "all-features run is CPU test evidence from this host only, not hosted CI nor GPU runtime evidence",\n'
  printf '  "started": "%s",\n' "$START_TS"
  printf '  "finished": "%s",\n' "$END_TS"
  printf '  "overall": "%s",\n' "$([ $OVERALL -eq 0 ] && echo pass || echo fail)"
  printf '  "gates": [\n'
  for i in "${!RESULTS[@]}"; do
    [[ $i -gt 0 ]] && printf ',\n'
    printf '    %s' "${RESULTS[$i]}"
  done
  printf '\n  ]\n'
  printf '}\n'
} > "$EVIDENCE_FILE"

echo ""
echo "Evidence written to $EVIDENCE_FILE"
echo "Overall: $([ $OVERALL -eq 0 ] && echo PASS || echo FAIL)"
exit $OVERALL
