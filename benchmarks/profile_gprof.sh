#!/usr/bin/env bash
# Profile the existing SYNSPEC 54 source with GNU gprof.
# Usage: bash benchmarks/profile_gprof.sh /path/to/input-case [shell command]
# Input directory must contain all files needed by the SYNSPEC run.
# Optional command runs inside copied case; default: './synspec54 < fort.5'
set -euo pipefail
if [[ $# -lt 1 ]]; then echo "Usage: $0 INPUT_CASE_DIR [COMMAND]" >&2; exit 2; fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE="$(cd "$1" && pwd)"
CMD="${2:-./synspec54 < fort.5}"
OUT="${PROFILE_OUTPUT:-$ROOT/results/gprof-$(date +%Y%m%d-%H%M%S)}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
if [[ "$OUT/" == "$CASE/"* || "$CASE/" == "$OUT/"* ]]; then
  echo "Profile output and input case must be separate directories" >&2; exit 2
fi
mkdir -p "$OUT/build"
echo "Compiler: $(gfortran --version | head -1)" > "$OUT/metadata.txt"
echo "Command: $CMD" >> "$OUT/metadata.txt"
echo "Source SHA256: $(sha256sum "$ROOT/src/synspec54.f")" >> "$OUT/metadata.txt"
flags=(-fno-automatic -O2 -g -pg -fno-omit-frame-pointer)
if [[ "$(uname -s)" == Linux ]]; then flags+=(-mcmodel=medium); fi
printf 'Flags: %s\n' "${flags[*]}" >> "$OUT/metadata.txt"
(cd "$ROOT/src" && gfortran "${flags[@]}" synspec54.f -o "$OUT/build/synspec54") \
  > "$OUT/build/compile.stdout" 2> "$OUT/build/compile.stderr"
mkdir -p "$OUT/run"
cp -a "$CASE/." "$OUT/run/"
cp "$OUT/build/synspec54" "$OUT/run/synspec54"
(
 cd "$OUT/run"
 export OMP_NUM_THREADS=1
 /usr/bin/time -p sh -c "$CMD"
) > "$OUT/run.stdout" 2> "$OUT/run.stderr"
if [[ ! -f "$OUT/run/gmon.out" ]]; then
 echo "No gmon.out generated: verify that SYNSPEC exited normally." >&2; exit 1
fi
gprof "$OUT/build/synspec54" "$OUT/run/gmon.out" > "$OUT/gprof.txt"
echo "Profile report: $OUT/gprof.txt"
echo "Inspect the flat profile for exclusive time and call graph for inclusive time."
