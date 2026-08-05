#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="${PIC_APP_DIR:-$repo_root/PIC-IFE_GEC}"
archive_root="${PIC_ARCHIVE_ROOT:-$repo_root/verification_runs}"
config_file="$app_dir/case_config.txt"

config_value() {
  local key="$1"
  awk -F= -v key="$key" '
    $1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
      value=$2
      sub(/^[[:space:]]+/, "", value)
      sub(/[[:space:]]+$/, "", value)
      print value
      exit
    }
  ' "$config_file"
}

if [[ ! -f "$config_file" ]]; then
  echo "missing case configuration: $config_file" >&2
  exit 1
fi

case_name="$(config_value case_name)"
t30_step="$(config_value target_omega_pi_t_30_step)"
t50_step="$(config_value target_omega_pi_t_50_step)"

if [[ -z "$case_name" || ! "$t30_step" =~ ^[0-9]+$ || ! "$t50_step" =~ ^[0-9]+$ ]]; then
  echo "invalid case_name or target steps in $config_file" >&2
  exit 1
fi

printf -v t30_file_step '%06d' "$t30_step"
printf -v t50_file_step '%06d' "$t50_step"

archive_dir="$archive_root/$case_name"
if [[ -e "$archive_dir" ]]; then
  echo "archive already exists; refusing to mix results: $archive_dir" >&2
  exit 2
fi

required_files=(
  "$app_dir/run.log"
  "$app_dir/OUTPUT/global_diagnostics.csv"
  "$app_dir/OUTPUT/physics_parameter.inp"
  "$app_dir/OUTPUT/normalize.inp"
  "$app_dir/OUTPUT/Field/Average_x_${t30_file_step}.dat"
  "$app_dir/OUTPUT/Field/Average_x_${t50_file_step}.dat"
  "$app_dir/OUTPUT/Velocity/velocity_IJ_3${t30_file_step}.dat"
  "$app_dir/OUTPUT/Velocity/velocity_IJ_3${t50_file_step}.dat"
)

for file in "${required_files[@]}"; do
  if [[ ! -f "$file" ]]; then
    echo "required result is missing: $file" >&2
    exit 1
  fi
done

mkdir -p "$archive_dir"
cp "$config_file" "$app_dir/INPUT/pic.inp" "$app_dir/MCC_jw/input/controlflow.txt" "$archive_dir/"
cp "$app_dir/run.log" "$app_dir/OUTPUT/global_diagnostics.csv" \
   "$app_dir/OUTPUT/physics_parameter.inp" "$app_dir/OUTPUT/normalize.inp" "$archive_dir/"
cp "$app_dir/OUTPUT/Field/Average_x_${t30_file_step}.dat" \
   "$app_dir/OUTPUT/Field/Average_x_${t50_file_step}.dat" "$archive_dir/"
cp "$app_dir/OUTPUT/Velocity/velocity_IJ_3${t30_file_step}.dat" \
   "$app_dir/OUTPUT/Velocity/velocity_IJ_3${t50_file_step}.dat" "$archive_dir/"

if [[ -f "$app_dir/run_metadata.txt" ]]; then
  cp "$app_dir/run_metadata.txt" "$archive_dir/"
fi

echo "archived verification case: $archive_dir"
