#!/usr/bin/env python3
"""Draw PNG figures from a saved PIC verification case.

This script reads a case directory under ``verification_runs``. It does not
depend on ``PIC-IFE_GEC/OUTPUT`` after the run has been archived.

Expected flat case directory:
  verification_runs/<case_name>/
    Average_x_012000.dat
    Average_x_020000.dat
    velocity_IJ_3012000.dat
    velocity_IJ_3020000.dat
    physics_parameter.inp
    case_config.txt

Usage:
  cd ~/pic-
  python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 12000
  python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 20000

If the step is omitted, the latest available Average_x/velocity_IJ pair is used.
PNG files are written to:
  verification_runs/<case_name>/postprocessed/
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - user environment check
    raise SystemExit(
        "matplotlib is required. Try: python3 -m pip install --user matplotlib"
    ) from exc


FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")
FIELD_RE = re.compile(r"Average_x_(\d{6})\.dat$")
VELOCITY_RE = re.compile(r"velocity_IJ_3(\d{6})\.dat$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw density and electron thermal-velocity PNGs from a saved verification case."
    )
    parser.add_argument(
        "case_dir",
        type=Path,
        help="Saved case directory, for example verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000",
    )
    parser.add_argument(
        "step",
        type=int,
        nargs="?",
        help="Iteration step such as 12000 or 20000. If omitted, use the latest output.",
    )
    parser.add_argument(
        "label",
        nargs="?",
        help="Optional title label, for example 'thermal reservoir'.",
    )
    return parser.parse_args()


def read_numeric_table(path: Path, min_cols: int) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            values: list[float] = []
            for token in line.replace(",", " ").split():
                try:
                    values.append(float(token.replace("D", "E").replace("d", "e")))
                except ValueError:
                    values = []
                    break
            if len(values) >= min_cols:
                rows.append(values)
    if not rows:
        raise SystemExit(f"No numeric rows found in {path}")
    return np.asarray(rows, dtype=float)


def read_velocity_table(path: Path) -> np.ndarray:
    """Read the seven-column velocity format, including wrapped records."""
    values: list[float] = []
    with path.open("r", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                row_values = [
                    float(token.replace("D", "E").replace("d", "e"))
                    for token in line.replace(",", " ").split()
                ]
            except ValueError:
                continue
            values.extend(row_values)

    column_count = 7
    if not values:
        raise SystemExit(f"No numeric velocity data found in {path}")
    if len(values) % column_count:
        raise SystemExit(
            f"Incomplete velocity record in {path}: "
            f"found {len(values)} values, expected a multiple of {column_count}"
        )
    return np.asarray(values, dtype=float).reshape((-1, column_count))


def parse_ref_value(path: Path, name: str, default: float) -> float:
    if not path.exists():
        return default
    for line in path.read_text(errors="replace").splitlines():
        if name.lower() in line.lower():
            nums = FLOAT_RE.findall(line)
            if nums:
                return float(nums[-1].replace("D", "E").replace("d", "e"))
    return default


def parse_config_value(path: Path, key: str, default: str | None = None) -> str | None:
    if not path.exists():
        return default
    prefix = key.lower()
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip().lower() == prefix:
            return rhs.strip()
    return default


def parse_config_float(path: Path, key: str, default: float) -> float:
    value = parse_config_value(path, key)
    if value is None:
        return default
    nums = FLOAT_RE.findall(value)
    if not nums:
        return default
    return float(nums[-1].replace("D", "E").replace("d", "e"))


def output_paths(case_dir: Path, pattern: str, nested: tuple[str, str]) -> list[Path]:
    paths = sorted(case_dir.glob(pattern))
    if paths:
        return paths
    subdir, subpattern = nested
    return sorted((case_dir / subdir).glob(subpattern))


def parse_step(path: Path, pattern: re.Pattern[str]) -> int | None:
    match = pattern.search(path.name)
    if not match:
        return None
    return int(match.group(1))


def choose_file(paths: list[Path], pattern: re.Pattern[str], requested: int | None) -> tuple[Path, int]:
    if not paths:
        raise SystemExit("No matching output files found in the case directory.")
    indexed: list[tuple[int, Path]] = []
    for path in paths:
        step = parse_step(path, pattern)
        if step is not None:
            indexed.append((step, path))
    if not indexed:
        raise SystemExit("Could not parse output step numbers.")
    indexed.sort()
    if requested is None:
        return indexed[-1][1], indexed[-1][0]
    exact = [item for item in indexed if item[0] == requested]
    if exact:
        return exact[0][1], exact[0][0]
    step, path = min(indexed, key=lambda item: abs(item[0] - requested))
    print(f"Requested step {requested:06d} not found; using nearest step {step:06d}: {path}")
    return path, step


def infer_label(config_file: Path, explicit_label: str | None) -> str:
    if explicit_label:
        return explicit_label
    left_boundary = parse_config_value(config_file, "left_boundary", "unknown")
    if left_boundary == "thermal":
        return "thermal reservoir"
    if left_boundary == "specular":
        return "specular reflection"
    return str(left_boundary)


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir
    if not case_dir.exists():
        raise SystemExit(f"Case directory does not exist: {case_dir}")

    field_file, field_step = choose_file(
        output_paths(case_dir, "Average_x_*.dat", ("OUTPUT/Field", "Average_x_*.dat")),
        FIELD_RE,
        args.step,
    )
    velocity_file, _ = choose_file(
        output_paths(case_dir, "velocity_IJ_3*.dat", ("OUTPUT/Velocity", "velocity_IJ_3*.dat")),
        VELOCITY_RE,
        field_step,
    )

    field = read_numeric_table(field_file, 9)
    velocity = read_velocity_table(velocity_file)

    physics_file = case_dir / "physics_parameter.inp"
    config_file = case_dir / "case_config.txt"
    if not physics_file.exists():
        physics_file = case_dir / "OUTPUT" / "physics_parameter.inp"

    n0 = parse_ref_value(physics_file, "density ref", parse_config_float(config_file, "density_ref_m3", 1.0e21))
    lambda_d0 = parse_ref_value(physics_file, "length ref", 1.0)
    mi_me = parse_config_float(config_file, "ion_mass_over_electron_mass", 400.0)
    dt_omega_pe = parse_config_float(config_file, "dt_omega_pe", 0.05)

    if lambda_d0 <= 0 or (lambda_d0 == 1.0 and np.nanmax(field[:, 0]) < 1.0):
        dx_values = np.diff(field[:, 0])
        dx_values = dx_values[np.isfinite(dx_values) & (dx_values > 0)]
        if dx_values.size:
            lambda_d0 = float(np.nanmedian(dx_values))

    x_lambda = field[:, 0] / lambda_d0
    ne = np.abs(field[:, 3]) / n0
    ni = np.abs(field[:, 4]) / n0

    dx_lambda_d = parse_config_float(config_file, "dx_lambdaD", 1.0)
    x_vel = (velocity[:, 6] - 0.5) * dx_lambda_d
    thermal_e = velocity[:, 2]
    count_e = velocity[:, 4]
    thermal_e = np.where(count_e > 0, thermal_e, np.nan)
    order = np.argsort(x_vel)
    x_vel = x_vel[order]
    thermal_e = thermal_e[order]

    front_mask = ni > 1.0e-3
    front_x = float(np.nanmax(x_lambda[front_mask])) if np.any(front_mask) else math.nan
    omega_pi_t = field_step * dt_omega_pe / math.sqrt(mi_me)

    out_dir = case_dir / "postprocessed"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"maxwell_mi400_t{field_step:06d}.png"

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=False)

    ax = axes[0]
    ax.plot(x_lambda, ne, label="electron", lw=1.8)
    ax.plot(x_lambda, ni, label="ion", lw=1.8)
    x_plot_max = min(350.0, float(np.nanmax(x_lambda)))
    if math.isfinite(front_x) and front_x <= x_plot_max:
        ax.axvline(front_x, color="0.35", lw=1.0, ls="--")
        ax.text(front_x, 0.04, f"front {front_x:.1f}", rotation=90, va="bottom", ha="right", fontsize=8)
    ax.set_xlim(0, x_plot_max)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("x / lambda_D0")
    ax.set_ylabel("n / n0")
    ax.set_title("Density")
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(x_vel, thermal_e, color="tab:red", lw=1.8)
    ax.set_xlim(0, min(350.0, float(np.nanmax(x_vel))))
    visible = thermal_e[(x_vel >= 0) & (x_vel <= 350.0) & np.isfinite(thermal_e)]
    if visible.size:
        ax.set_ylim(0, max(1.0, float(np.nanpercentile(visible, 99.0)) * 1.15))
    else:
        ax.set_ylim(bottom=0)
    ax.set_xlabel("x / lambda_D0")
    ax.set_ylabel("v_th,e / v_th,e0")
    ax.set_title("Electron thermal velocity")

    label = infer_label(config_file, args.label)
    fig.suptitle(
        f"Maxwellian, {label}, mi/me={mi_me:g}, step={field_step}, omega_pi t={omega_pi_t:.1f}",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.82, wspace=0.24)
    fig.savefig(out_file, dpi=220)
    print(f"field:     {field_file}")
    print(f"velocity:  {velocity_file}")
    print(f"figure:    {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
