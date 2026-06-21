#!/usr/bin/env python3
"""Plot the Maxwellian mi/me=400 plasma-expansion proof case.

Usage:
  python3 scripts/plot_maxwell_mi400.py [PIC-IFE_GEC] [step]

If step is omitted, the latest available Average_x/velocity_IJ output is used.
The figure is written to PIC-IFE_GEC/figures/.
"""

import math
import re
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "matplotlib is required for plotting. Try: python3 -m pip install --user matplotlib"
    ) from exc


FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")


def read_numeric_table(path, min_cols):
    rows = []
    with path.open("r", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            values = []
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


def parse_ref_value(path, name, default):
    if not path.exists():
        return default
    for line in path.read_text(errors="replace").splitlines():
        if name.lower() in line.lower():
            nums = FLOAT_RE.findall(line)
            if nums:
                return float(nums[-1].replace("D", "E").replace("d", "e"))
    return default


def parse_step(path, pattern):
    match = pattern.search(path.name)
    if not match:
        return None
    return int(match.group(1))


def choose_file(paths, pattern, requested):
    if not paths:
        raise SystemExit("No matching output files found. Did the simulation finish and write diagnostics?")
    indexed = []
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


def main() -> int:
    case_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("PIC-IFE_GEC")
    requested_step = int(sys.argv[2]) if len(sys.argv) > 2 else None
    boundary_label = sys.argv[3] if len(sys.argv) > 3 else "specular reflection"

    field_pattern = re.compile(r"Average_x_(\d{6})\.dat$")
    velocity_pattern = re.compile(r"velocity_IJ_3(\d{6})\.dat$")

    field_file, field_step = choose_file(
        sorted((case_dir / "OUTPUT" / "Field").glob("Average_x_*.dat")),
        field_pattern,
        requested_step,
    )
    velocity_file, _ = choose_file(
        sorted((case_dir / "OUTPUT" / "Velocity").glob("velocity_IJ_3*.dat")),
        velocity_pattern,
        field_step,
    )

    field = read_numeric_table(field_file, 9)
    velocity = read_numeric_table(velocity_file, 7)

    physics_file = case_dir / "OUTPUT" / "physics_parameter.inp"
    lambda_d0 = parse_ref_value(physics_file, "length ref", 1.0)
    n0 = parse_ref_value(physics_file, "density ref", 1.0e21)

    if lambda_d0 <= 0 or (lambda_d0 == 1.0 and np.nanmax(field[:, 0]) < 1.0):
        dx_values = np.diff(field[:, 0])
        dx_values = dx_values[np.isfinite(dx_values) & (dx_values > 0)]
        if dx_values.size:
            lambda_d0 = float(np.nanmedian(dx_values))

    x_lambda = field[:, 0] / lambda_d0
    ne = np.abs(field[:, 3]) / n0
    ni = np.abs(field[:, 4]) / n0

    if velocity.shape[1] >= 10:
        x_vel = velocity[:, 9] - 0.5
        thermal_e = velocity[:, 3]
        count_e = velocity[:, 6]
    else:
        x_vel = velocity[:, 6] - 0.5
        thermal_e = velocity[:, 2]
        count_e = velocity[:, 4]
    thermal_e = np.where(count_e > 0, thermal_e, np.nan)
    order = np.argsort(x_vel)
    x_vel = x_vel[order]
    thermal_e = thermal_e[order]

    front_mask = ni > 1.0e-3
    front_x = float(np.nanmax(x_lambda[front_mask])) if np.any(front_mask) else math.nan
    omega_pi_t = field_step * 0.05 / math.sqrt(400.0)

    out_dir = case_dir / "figures"
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

    fig.suptitle(f"Maxwellian, {boundary_label}, mi/me=400, step={field_step}, omega_pi t={omega_pi_t:.1f}", fontsize=13)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.82, wspace=0.24)
    fig.savefig(out_file, dpi=220)
    print(f"field:    {field_file}")
    print(f"velocity: {velocity_file}")
    print(f"figure:   {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
