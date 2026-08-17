#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/setup_paper_case.sh DISTRIBUTION SHAPE PPC NT BOUNDARY SEED DT DX MI_ME
#
# Examples:
#   bash scripts/setup_paper_case.sh maxwellian none 1000 20000 thermal 101 0.05 1.0 400
#   bash scripts/setup_paper_case.sh kappa 2 1000 20000 thermal 101 0.05 1.0 400
#   bash scripts/setup_paper_case.sh polytropic 2 1000 20000 thermal 101 0.05 1.0 400

if [[ "$#" -ne 9 ]]; then
  echo "usage: $0 DISTRIBUTION SHAPE PPC NT BOUNDARY SEED DT DX MI_ME" >&2
  exit 2
fi

distribution="$1"
shape="$2"
ppc="$3"
nt="$4"
boundary="$5"
seed="$6"
dt="$7"
dx="$8"
mass_ratio="$9"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$script_dir/setup_maxwell_mi400_case.sh" \
  "$ppc" "$nt" "$boundary" "$seed" "$dt" "$dx" \
  "$distribution" "$shape" "$mass_ratio"
