#!/usr/bin/env bash
# Offline fixture check for local-ci.sh SHA/dirty handling. No cargo is run:
# the gates fail fast against a stub PATH, we only test pre-gate refusals.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T" && git init -q . && git -c user.email=a@b -c user.name=t commit -q --allow-empty -m a
cp "$ROOT/local-ci.sh" .
mkdir -p "$T.bin"; printf "#!/bin/sh\nexit 0\n" > "$T.bin/cargo"; chmod +x "$T.bin/cargo"
export PATH="$T.bin:$PATH"; trap 'rm -rf "$T" "$T.bin"' EXIT
fail=0
chk() { [[ "$2" == "$3" ]] || { echo "FAIL $1: got $2 want $3"; fail=1; }; }
git -c user.email=a@b -c user.name=t commit -q --allow-empty -m b
./local-ci.sh --sha HEAD~1 >/dev/null 2>&1; chk "foreign sha refused" $? 2
./local-ci.sh --sha HEAD >/dev/null 2>&1; chk "untracked wrapper counts dirty" $? 2
git add local-ci.sh; git -c user.email=a@b -c user.name=t commit -q -m c
./local-ci.sh --sha "$(git rev-parse HEAD)" >/dev/null 2>&1; chk "clean sha ok" $? 0
grep -q '"sha_bound": true' local-ci-evidence.json; chk "bound evidence" $? 0
echo x > dirty.txt
./local-ci.sh >/dev/null 2>&1; chk "dirty refused" $? 2
mkdir -p subdir
(cd subdir && ../local-ci.sh >/dev/null 2>&1); chk "root dirty refused from subdirectory" $? 2
./local-ci.sh --allow-dirty >/dev/null 2>&1; chk "allow-dirty forces fail" $? 1
grep -q '"sha_bound": false' local-ci-evidence.json; chk "unbound evidence" $? 0
[[ $fail -eq 0 ]] && echo "local-ci wrapper checks passed"; exit $fail
