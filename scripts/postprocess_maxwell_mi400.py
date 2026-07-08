#!/usr/bin/env python3
"""Postprocess a saved Maxwellian mi/me=400 verification case.

This script reads archived files from ``verification_runs/<case_name>`` and
writes all derived products to ``verification_runs/<case_name>/postprocessed``.
It does not require ``PIC-IFE_GEC/OUTPUT`` after the run has been saved.

Usage:
  cd ~/pic-
  python3 scripts/postprocess_maxwell_mi400.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000

Optional custom steps:
  python3 scripts/postprocess_maxwell_mi400.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 12000 20000

Generated files:
  postprocessed/profiles_t30.csv
  postprocessed/profiles_t50.csv
  postprocessed/gamma_fit_t30.csv
  postprocessed/gamma_fit_t50.csv
  postprocessed/postprocess_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")
FIELD_RE = re.compile(r"Average_x_(\d{6})\.dat$")
VELOCITY_RE = re.compile(r"velocity_IJ_3(\d{6})\.dat$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create profile CSVs and gamma_e fits from a saved verification case."
    )
    parser.add_argument(
        "case_dir",
        type=Path,
        help="Saved case directory, for example verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000",
    )
    parser.add_argument(
        "steps",
        type=int,
        nargs="*",
        help="Optional steps such as 12000 20000. Defaults to t30/t50 from case_config.txt.",
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


def parse_config_int(path: Path, key: str, default: int) -> int:
    return int(round(parse_config_float(path, key, float(default))))


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


def choose_file(paths: list[Path], pattern: re.Pattern[str], requested: int) -> tuple[Path, int]:
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
    exact = [item for item in indexed if item[0] == requested]
    if exact:
        return exact[0][1], exact[0][0]
    step, path = min(indexed, key=lambda item: abs(item[0] - requested))
    print(f"Requested step {requested:06d} not found; using nearest step {step:06d}: {path}")
    return path, step


def label_for_step(config_file: Path, step: int, omega_pi_t: float) -> str:
    t30_step = parse_config_int(config_file, "target_omega_pi_t_30_step", 12000)
    t50_step = parse_config_int(config_file, "target_omega_pi_t_50_step", 20000)
    if step == t30_step:
        return "t30"
    if step == t50_step:
        return "t50"
    return f"step{step:06d}_opit{omega_pi_t:.1f}".replace(".", "p")


def build_profile(
    field_file: Path,
    velocity_file: Path,
    physics_file: Path,
    config_file: Path,
) -> tuple[np.ndarray, dict[str, float]]:
    field = read_numeric_table(field_file, 9)
    velocity = read_numeric_table(velocity_file, 10)

    n0 = parse_ref_value(physics_file, "density ref", parse_config_float(config_file, "density_ref_m3", 1.0e21))
    lambda_d0 = parse_ref_value(physics_file, "length ref", 1.0)
    if lambda_d0 <= 0 or (lambda_d0 == 1.0 and np.nanmax(field[:, 0]) < 1.0):
        dx_values = np.diff(field[:, 0])
        dx_values = dx_values[np.isfinite(dx_values) & (dx_values > 0)]
        if dx_values.size:
            lambda_d0 = float(np.nanmedian(dx_values))

    x_field = field[:, 0] / lambda_d0
    phi = field[:, 1]
    rho = field[:, 2] / n0
    ne = np.abs(field[:, 3]) / n0
    ni = np.abs(field[:, 4]) / n0
    ek_e = field[:, 5]
    ek_i = field[:, 6]

    x_cell = velocity[:, 9] - 0.5
    drift_e = velocity[:, 0]
    drift_i = velocity[:, 1]
    vthe = velocity[:, 3]
    vthi = velocity[:, 4]
    count_e = velocity[:, 6]
    count_i = velocity[:, 7]

    order = np.argsort(x_cell)
    x_cell = x_cell[order]
    drift_e = drift_e[order]
    drift_i = drift_i[order]
    vthe = vthe[order]
    vthi = vthi[order]
    count_e = count_e[order]
    count_i = count_i[order]

    ne_cell = np.interp(x_cell, x_field, ne)
    ni_cell = np.interp(x_cell, x_field, ni)
    rho_cell = np.interp(x_cell, x_field, rho)
    phi_cell = np.interp(x_cell, x_field, phi)
    ek_e_cell = np.interp(x_cell, x_field, ek_e)
    ek_i_cell = np.interp(x_cell, x_field, ek_i)

    vthe = np.where(count_e > 0, vthe, np.nan)
    vthi = np.where(count_i > 0, vthi, np.nan)
    te = vthe * vthe
    ti = vthi * vthi

    profile = np.column_stack(
        [
            x_cell,
            ne_cell,
            ni_cell,
            rho_cell,
            phi_cell,
            vthe,
            vthi,
            te,
            ti,
            drift_e,
            drift_i,
            count_e,
            count_i,
            ek_e_cell,
            ek_i_cell,
        ]
    )
    metadata = {
        "n0": n0,
        "lambda_d0": lambda_d0,
        "profile_rows": float(profile.shape[0]),
    }
    return profile, metadata


def fit_gamma(profile: np.ndarray) -> dict[str, float]:
    x = profile[:, 0]
    ne = profile[:, 1]
    te = profile[:, 7]
    count_e = profile[:, 11]

    mask = (
        np.isfinite(x)
        & np.isfinite(ne)
        & np.isfinite(te)
        & (ne > 1.0e-3)
        & (ne < 0.95)
        & (te > 0.0)
        & (count_e >= 10.0)
    )
    if int(np.count_nonzero(mask)) < 5:
        mask = (
            np.isfinite(ne)
            & np.isfinite(te)
            & (ne > 1.0e-6)
            & (ne < 1.0)
            & (te > 0.0)
            & (count_e > 0.0)
        )
    n_points = int(np.count_nonzero(mask))
    if n_points < 3:
        return {
            "gamma_e": math.nan,
            "slope": math.nan,
            "intercept": math.nan,
            "r2": math.nan,
            "n_points": float(n_points),
            "ne_min": math.nan,
            "ne_max": math.nan,
            "te_min": math.nan,
            "te_max": math.nan,
        }

    log_ne = np.log(ne[mask])
    log_te = np.log(te[mask])
    slope, intercept = np.polyfit(log_ne, log_te, 1)
    predicted = slope * log_ne + intercept
    ss_res = float(np.sum((log_te - predicted) ** 2))
    ss_tot = float(np.sum((log_te - np.mean(log_te)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else math.nan
    return {
        "gamma_e": float(slope + 1.0),
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": r2,
        "n_points": float(n_points),
        "ne_min": float(np.nanmin(ne[mask])),
        "ne_max": float(np.nanmax(ne[mask])),
        "te_min": float(np.nanmin(te[mask])),
        "te_max": float(np.nanmax(te[mask])),
    }


def write_profile(path: Path, profile: np.ndarray) -> None:
    header = [
        "x_over_lambda_D0",
        "ne_over_n0",
        "ni_over_n0",
        "rho_over_n0",
        "phi",
        "vth_e_over_vth_e0",
        "vth_i_over_vth_e0",
        "Te_over_Te0",
        "Ti_over_Te0",
        "drift_e",
        "drift_i",
        "electron_count",
        "ion_count",
        "Ek_ele",
        "Ek_ion",
    ]
    np.savetxt(path, profile, delimiter=",", header=",".join(header), comments="")


def write_fit(path: Path, row: dict[str, str | float | int]) -> None:
    fields = [
        "step",
        "omega_pi_t",
        "gamma_e",
        "slope",
        "intercept",
        "r2",
        "n_points",
        "ne_min",
        "ne_max",
        "te_min",
        "te_max",
        "field_file",
        "velocity_file",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir
    if not case_dir.exists():
        raise SystemExit(f"Case directory does not exist: {case_dir}")

    config_file = case_dir / "case_config.txt"
    physics_file = case_dir / "physics_parameter.inp"
    if not physics_file.exists():
        physics_file = case_dir / "OUTPUT" / "physics_parameter.inp"

    mi_me = parse_config_float(config_file, "ion_mass_over_electron_mass", 400.0)
    dt_omega_pe = parse_config_float(config_file, "dt_omega_pe", 0.05)
    steps = args.steps or [
        parse_config_int(config_file, "target_omega_pi_t_30_step", 12000),
        parse_config_int(config_file, "target_omega_pi_t_50_step", 20000),
    ]

    field_paths = output_paths(case_dir, "Average_x_*.dat", ("OUTPUT/Field", "Average_x_*.dat"))
    velocity_paths = output_paths(case_dir, "velocity_IJ_3*.dat", ("OUTPUT/Velocity", "velocity_IJ_3*.dat"))
    out_dir = case_dir / "postprocessed"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_lines = [
        f"case_dir = {case_dir}",
        f"postprocessed_dir = {out_dir}",
        f"mi_me = {mi_me:g}",
        f"dt_omega_pe = {dt_omega_pe:g}",
    ]

    for requested in steps:
        field_file, step = choose_file(field_paths, FIELD_RE, requested)
        velocity_file, velocity_step = choose_file(velocity_paths, VELOCITY_RE, step)
        omega_pi_t = step * dt_omega_pe / math.sqrt(mi_me)
        label = label_for_step(config_file, step, omega_pi_t)

        profile, metadata = build_profile(field_file, velocity_file, physics_file, config_file)
        fit = fit_gamma(profile)

        profile_file = out_dir / f"profiles_{label}.csv"
        fit_file = out_dir / f"gamma_fit_{label}.csv"
        write_profile(profile_file, profile)
        write_fit(
            fit_file,
            {
                "step": step,
                "omega_pi_t": omega_pi_t,
                "gamma_e": fit["gamma_e"],
                "slope": fit["slope"],
                "intercept": fit["intercept"],
                "r2": fit["r2"],
                "n_points": int(fit["n_points"]),
                "ne_min": fit["ne_min"],
                "ne_max": fit["ne_max"],
                "te_min": fit["te_min"],
                "te_max": fit["te_max"],
                "field_file": field_file.name,
                "velocity_file": velocity_file.name,
            },
        )

        summary_lines.extend(
            [
                "",
                f"[{label}]",
                f"requested_step = {requested}",
                f"field_step = {step}",
                f"velocity_step = {velocity_step}",
                f"omega_pi_t = {omega_pi_t:.8g}",
                f"profile = {profile_file.name}",
                f"gamma_fit = {fit_file.name}",
                f"gamma_e = {fit['gamma_e']:.8g}",
                f"fit_r2 = {fit['r2']:.8g}",
                f"fit_points = {int(fit['n_points'])}",
                f"n0 = {metadata['n0']:.8e}",
                f"lambda_d0 = {metadata['lambda_d0']:.8e}",
            ]
        )

        print(f"{label}: wrote {profile_file}")
        print(f"{label}: wrote {fit_file}")

    summary_file = out_dir / "postprocess_summary.txt"
    summary_file.write_text("\n".join(summary_lines) + "\n")
    print(f"summary: wrote {summary_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
