#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$repo_root/build-velocity-distribution-test"

if [[ -n "${FC:-}" ]]; then
  compiler="$FC"
elif command -v ifort >/dev/null 2>&1; then
  compiler="ifort"
elif command -v ifx >/dev/null 2>&1; then
  compiler="ifx"
elif command -v gfortran >/dev/null 2>&1; then
  compiler="gfortran"
else
  echo "no Fortran compiler found; load intel/2022.1 or a GNU compiler module" >&2
  exit 1
fi

rm -rf "$build_dir"
mkdir -p "$build_dir"
cd "$build_dir"

"$compiler" -O2 -o velocity_distribution_smoke \
  "$repo_root/PIC-IFE_GEC/code/PIC/DRandomSeeded.f90" \
  "$repo_root/PIC-IFE_GEC/MCC_jw/code/base/VelocityDistribution.f90" \
  "$repo_root/tests/velocity_distribution_smoke.f90"

./velocity_distribution_smoke
