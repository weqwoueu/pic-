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
  postprocessed/gamma_fit_t30.png
  postprocessed/gamma_fit_t50.png
  postprocessed/postprocess_summary.txt
"""

from __future__ import annotations

import argparse
import csv
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
        "matplotlib is required. Try: python -m pip install --user matplotlib"
    ) from exc


FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")
FIELD_RE = re.compile(r"Average_x_(\d{6})\.dat$")
VELOCITY_RE = re.compile(r"velocity_IJ_3(\d{6})\.dat$")
FIT_NI_MIN = 1.0e-3
FIT_NI_MAX = 1.0
FIT_MIN_ELECTRON_COUNT = 10.0


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
    velocity = read_velocity_table(velocity_file)

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

    dx_lambda_d = parse_config_float(config_file, "dx_lambdaD", 1.0)
    x_cell = (velocity[:, 6] - 0.5) * dx_lambda_d
    drift_e = velocity[:, 0]
    drift_i = velocity[:, 1]
    vthe = velocity[:, 2]
    vthi = velocity[:, 3]
    count_e = velocity[:, 4]
    count_i = velocity[:, 5]

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


def fit_gamma(profile: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
    x = profile[:, 0]
    ne = profile[:, 1]
    ni = profile[:, 2]
    te = profile[:, 7]
    count_e = profile[:, 11]

    # The advisor's fixed fit interval is defined by ion density, while the
    # fitted thermodynamic relation remains Te(ne).
    mask = (
        np.isfinite(x)
        & np.isfinite(ne)
        & np.isfinite(ni)
        & np.isfinite(te)
        & (ni >= FIT_NI_MIN)
        & (ni < FIT_NI_MAX)
        & (ne > 0.0)
        & (te > 0.0)
        & (count_e >= FIT_MIN_ELECTRON_COUNT)
    )
    n_points = int(np.count_nonzero(mask))
    if n_points < 3:
        return {
            "gamma_e": math.nan,
            "gamma_e_standard_error": math.nan,
            "slope": math.nan,
            "slope_standard_error": math.nan,
            "intercept": math.nan,
            "intercept_standard_error": math.nan,
            "r2": math.nan,
            "residual_standard_error": math.nan,
            "n_points": float(n_points),
            "excluded_points": float(profile.shape[0] - n_points),
            "x_min": math.nan,
            "x_max": math.nan,
            "ni_min": math.nan,
            "ni_max": math.nan,
            "ne_min": math.nan,
            "ne_max": math.nan,
            "te_min": math.nan,
            "te_max": math.nan,
        }, mask

    log_ne = np.log(ne[mask])
    log_te = np.log(te[mask])
    slope, intercept = np.polyfit(log_ne, log_te, 1)
    predicted = slope * log_ne + intercept
    ss_res = float(np.sum((log_te - predicted) ** 2))
    ss_tot = float(np.sum((log_te - np.mean(log_te)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else math.nan
    degrees_of_freedom = n_points - 2
    residual_variance = ss_res / degrees_of_freedom
    residual_standard_error = math.sqrt(residual_variance)
    centered_log_ne = log_ne - np.mean(log_ne)
    sxx = float(np.sum(centered_log_ne**2))
    if sxx > 0.0:
        slope_standard_error = math.sqrt(residual_variance / sxx)
        intercept_standard_error = math.sqrt(
            residual_variance * (1.0 / n_points + float(np.mean(log_ne)) ** 2 / sxx)
        )
    else:
        slope_standard_error = math.nan
        intercept_standard_error = math.nan

    return {
        "gamma_e": float(slope + 1.0),
        "gamma_e_standard_error": slope_standard_error,
        "slope": float(slope),
        "slope_standard_error": slope_standard_error,
        "intercept": float(intercept),
        "intercept_standard_error": intercept_standard_error,
        "r2": r2,
        "residual_standard_error": residual_standard_error,
        "n_points": float(n_points),
        "excluded_points": float(profile.shape[0] - n_points),
        "x_min": float(np.nanmin(x[mask])),
        "x_max": float(np.nanmax(x[mask])),
        "ni_min": float(np.nanmin(ni[mask])),
        "ni_max": float(np.nanmax(ni[mask])),
        "ne_min": float(np.nanmin(ne[mask])),
        "ne_max": float(np.nanmax(ne[mask])),
        "te_min": float(np.nanmin(te[mask])),
        "te_max": float(np.nanmax(te[mask])),
    }, mask


def write_profile(path: Path, profile: np.ndarray, fit_mask: np.ndarray) -> None:
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
        "gamma_fit_included",
    ]
    output = np.column_stack([profile, fit_mask.astype(int)])
    np.savetxt(path, output, delimiter=",", header=",".join(header), comments="")


def write_fit(path: Path, row: dict[str, str | float | int]) -> None:
    fields = [
        "step",
        "omega_pi_t",
        "gamma_e",
        "gamma_e_standard_error",
        "slope",
        "slope_standard_error",
        "intercept",
        "intercept_standard_error",
        "r2",
        "residual_standard_error",
        "n_points",
        "excluded_points",
        "fit_ni_lower_inclusive",
        "fit_ni_upper_exclusive",
        "minimum_electron_count",
        "fit_log_base",
        "x_min",
        "x_max",
        "ni_min",
        "ni_max",
        "ne_min",
        "ne_max",
        "te_min",
        "te_max",
        "tail_exclusion",
        "field_file",
        "velocity_file",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)


def write_fit_figure(
    path: Path,
    profile: np.ndarray,
    fit_mask: np.ndarray,
    fit: dict[str, float],
    omega_pi_t: float,
) -> None:
    ne = profile[:, 1]
    te = profile[:, 7]
    count_e = profile[:, 11]
    valid = (
        np.isfinite(ne)
        & np.isfinite(te)
        & (ne > 0.0)
        & (te > 0.0)
        & (count_e >= FIT_MIN_ELECTRON_COUNT)
    )
    excluded = valid & ~fit_mask

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    if np.any(excluded):
        ax.scatter(
            np.log10(ne[excluded]),
            np.log10(te[excluded]),
            s=8,
            color="0.75",
            alpha=0.45,
            label="excluded",
        )
    if np.any(fit_mask):
        fit_x = np.log10(ne[fit_mask])
        fit_y = np.log10(te[fit_mask])
        ax.scatter(fit_x, fit_y, s=12, color="tab:blue", alpha=0.75, label="fit points")
        line_x = np.linspace(float(np.min(fit_x)), float(np.max(fit_x)), 200)
        line_y = fit["slope"] * line_x + fit["intercept"] / math.log(10.0)
        ax.plot(line_x, line_y, color="tab:red", lw=1.8, label="linear fit")

    ax.set_xlabel(r"$\log_{10}(n_e/n_{e0})$")
    ax.set_ylabel(r"$\log_{10}(T_e/T_{e0})$")
    ax.set_title(rf"$\omega_{{pi}}t={omega_pi_t:g}$")
    ax.text(
        0.04,
        0.96,
        rf"$\gamma_e={fit['gamma_e']:.4f}\pm{fit['gamma_e_standard_error']:.4f}$" "\n"
        rf"$R^2={fit['r2']:.4f}$, N = {int(fit['n_points'])}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)


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
        fit, fit_mask = fit_gamma(profile)

        profile_file = out_dir / f"profiles_{label}.csv"
        fit_file = out_dir / f"gamma_fit_{label}.csv"
        fit_figure = out_dir / f"gamma_fit_{label}.png"
        write_profile(profile_file, profile, fit_mask)
        write_fit(
            fit_file,
            {
                "step": step,
                "omega_pi_t": omega_pi_t,
                "gamma_e": fit["gamma_e"],
                "gamma_e_standard_error": fit["gamma_e_standard_error"],
                "slope": fit["slope"],
                "slope_standard_error": fit["slope_standard_error"],
                "intercept": fit["intercept"],
                "intercept_standard_error": fit["intercept_standard_error"],
                "r2": fit["r2"],
                "residual_standard_error": fit["residual_standard_error"],
                "n_points": int(fit["n_points"]),
                "excluded_points": int(fit["excluded_points"]),
                "fit_ni_lower_inclusive": FIT_NI_MIN,
                "fit_ni_upper_exclusive": FIT_NI_MAX,
                "minimum_electron_count": FIT_MIN_ELECTRON_COUNT,
                "fit_log_base": "natural",
                "x_min": fit["x_min"],
                "x_max": fit["x_max"],
                "ni_min": fit["ni_min"],
                "ni_max": fit["ni_max"],
                "ne_min": fit["ne_min"],
                "ne_max": fit["ne_max"],
                "te_min": fit["te_min"],
                "te_max": fit["te_max"],
                "tail_exclusion": "none",
                "field_file": field_file.name,
                "velocity_file": velocity_file.name,
            },
        )
        write_fit_figure(fit_figure, profile, fit_mask, fit, omega_pi_t)

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
                f"gamma_fit_figure = {fit_figure.name}",
                "fit_relation = ln(Te/Te0) = (gamma_e - 1) ln(ne/ne0) + b",
                f"fit_interval = {FIT_NI_MIN:g} <= ni/n0 < {FIT_NI_MAX:g}",
                f"minimum_electron_count = {FIT_MIN_ELECTRON_COUNT:g}",
                "tail_exclusion = none",
                f"gamma_e = {fit['gamma_e']:.8g}",
                f"gamma_e_standard_error = {fit['gamma_e_standard_error']:.8g}",
                f"fit_r2 = {fit['r2']:.8g}",
                f"residual_standard_error = {fit['residual_standard_error']:.8g}",
                f"fit_points = {int(fit['n_points'])}",
                f"excluded_points = {int(fit['excluded_points'])}",
                f"fit_x_range = [{fit['x_min']:.8g}, {fit['x_max']:.8g}]",
                f"n0 = {metadata['n0']:.8e}",
                f"lambda_d0 = {metadata['lambda_d0']:.8e}",
            ]
        )

        print(f"{label}: wrote {profile_file}")
        print(f"{label}: wrote {fit_file}")
        print(f"{label}: wrote {fit_figure}")

    summary_file = out_dir / "postprocess_summary.txt"
    summary_file.write_text("\n".join(summary_lines) + "\n")
    print(f"summary: wrote {summary_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
