#!/usr/bin/env bash
# local-ci.sh — Replay the full CI suite locally, including all-features tests.
#
# Usage:
#   ./local-ci.sh [--sha <SHA>]   # defaults to HEAD
#
# Requirements: Rust stable toolchain with rustfmt + clippy components.
#   rustup component add rustfmt clippy
#
# Exit code: 0 = all gates passed, non-zero = at least one gate failed.
#
# Evidence record is written to local-ci-evidence.json at the repo root.
# Each run is SHA-bound; the file is gitignored.

set -euo pipefail

SHA="$(git rev-parse HEAD 2>/dev/null || echo "unknown")"
case "$#" in
  0) ;;
  2)
    if [[ "$1" != "--sha" || -z "$2" ]]; then
      echo "Usage: $0 [--sha <SHA>]" >&2
      exit 2
    fi
    SHA="$2"
    ;;
  *)
    echo "Usage: $0 [--sha <SHA>]" >&2
    exit 2
    ;;
esac

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
echo "SHA: $SHA"
echo "Started: $START_TS"
echo ""

OVERALL=0
run_gate "rustfmt"          cargo fmt --all -- --check          || OVERALL=1
run_gate "clippy"           cargo clippy --all-targets -- -D warnings || OVERALL=1
run_gate "test-default"     cargo test                          || OVERALL=1
run_gate "test-all-features" cargo test --all-features          || OVERALL=1

END_TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

# Write evidence JSON
{
  printf '{\n'
  printf '  "schema": 1,\n'
  printf '  "sha": "%s",\n' "$SHA"
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
