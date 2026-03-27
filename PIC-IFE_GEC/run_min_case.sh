#!/usr/bin/env bash
set -euo pipefail

# One-command minimal case runner for Xibei Cluster (or any Linux with cmake+ifort/ifx).
# Defaults target a quick validation run.
NT="${NT:-1000}"
DT="${DT:-0.05}"

if [[ ! -f "CMakeLists.txt" || ! -f "INPUT/pic.inp" ]]; then
  echo "ERROR: run this script inside PIC-IFE_GEC/." >&2
  exit 1
fi

cp -f INPUT/pic.inp INPUT/pic.inp.bak

python - <<'PY'
from pathlib import Path
import os, re, sys

p = Path("INPUT/pic.inp")
s = p.read_text(encoding="utf-8", errors="ignore")
nt = os.environ.get("NT", "1000")
dt = os.environ.get("DT", "0.05")

# Replace first active "nt, dt" line only.
pat = re.compile(r"(?m)^(\s*)(\d+)\s*,\s*([0-9.+\-DdEe]+)(\s*!.*nt\s*,\s*dt.*)?\s*$")
m = pat.search(s)
if not m:
    print("ERROR: cannot find a valid 'nt, dt' line in INPUT/pic.inp", file=sys.stderr)
    sys.exit(2)

indent = m.group(1)
comment = m.group(4) or " !nt, dt"
new_line = f"{indent}{nt}, {dt}{comment}"
s2 = s[:m.start()] + new_line + s[m.end():]
p.write_text(s2, encoding="utf-8")
print(f"Updated nt,dt -> {new_line}")
PY

mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP
rm -rf build

FC="${FC:-ifort}"
if ! command -v "$FC" >/dev/null 2>&1; then
  if command -v ifx >/dev/null 2>&1; then
    FC="ifx"
  else
    echo "ERROR: ifort/ifx not found. Load Intel compiler first." >&2
    exit 3
  fi
fi

echo "Using Fortran compiler: $FC"
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(command -v "$FC")"
cmake --build build -j"$(nproc)"
test -x ./1DPIC
./1DPIC
