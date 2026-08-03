#!/usr/bin/env python3
"""Summarize gamma_e results from repeated Maxwellian PIC seeds.

Run this after postprocess_maxwell_mi400.py has been applied to every case.

Example:
  python scripts/summarize_maxwell_seed_ensemble.py \
    verification_runs/maxwellian_dx1_dt005_ppc80000_seed101_thermal \
    verification_runs/maxwellian_dx1_dt005_ppc80000_seed202_thermal \
    verification_runs/maxwellian_dx1_dt005_ppc80000_seed303_thermal \
    --reference 1.023 --reference-error 0.003

Generated files are placed in a derived seed-ensemble directory unless
--output-dir is supplied:
  postprocessed/seed_results.csv
  postprocessed/seed_statistics.csv
  postprocessed/seed_statistics.txt
  postprocessed/seed_gamma.png
"""

import argparse
import csv
import math
import re
import statistics
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - environment check
    raise SystemExit(
        "matplotlib is required. Try: python -m pip install --user matplotlib"
    ) from exc


SEED_RE = re.compile(r"_seed(\d+)_")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute equal-weight seed statistics from postprocessed cases."
    )
    parser.add_argument(
        "case_dirs",
        type=Path,
        nargs="+",
        help="Saved verification case directories with postprocessed summaries.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Optional output directory. Defaults to a sibling seed-ensemble case.",
    )
    parser.add_argument(
        "--reference",
        type=float,
        help="Optional literature/reference gamma_e value to show in the figure.",
    )
    parser.add_argument(
        "--reference-error",
        type=float,
        help="Optional symmetric uncertainty for --reference.",
    )
    return parser.parse_args()


def parse_key_value_summary(path):
    general = {}
    sections = {}
    current = general
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                label = line[1:-1]
                current = sections.setdefault(label, {})
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                current[key.strip()] = value.strip()
    return general, sections


def seed_from_case(case_dir):
    match = SEED_RE.search(case_dir.name)
    if not match:
        raise SystemExit(f"Cannot find _seedNNN_ in case name: {case_dir.name}")
    return int(match.group(1))


def canonical_case_name(case_dir):
    return SEED_RE.sub("_seed_ensemble_", case_dir.name)


def default_output_dir(case_dirs):
    name = canonical_case_name(case_dirs[0])
    return case_dirs[0].parent / name / "postprocessed"


def require_float(mapping, key, source):
    try:
        return float(mapping[key])
    except (KeyError, ValueError) as exc:
        raise SystemExit(f"Missing or invalid {key} in {source}") from exc


def load_cases(case_dirs):
    records = []
    expected_name = canonical_case_name(case_dirs[0])
    expected_general = None
    expected_labels = None
    seen_seeds = set()

    for case_dir in case_dirs:
        if canonical_case_name(case_dir) != expected_name:
            raise SystemExit("All case names must differ only by the seed number")
        seed = seed_from_case(case_dir)
        if seed in seen_seeds:
            raise SystemExit(f"Duplicate seed: {seed}")
        seen_seeds.add(seed)

        summary_path = case_dir / "postprocessed" / "postprocess_summary.txt"
        if not summary_path.is_file():
            raise SystemExit(
                f"Missing {summary_path}; run postprocess_maxwell_mi400.py first"
            )
        general, sections = parse_key_value_summary(summary_path)
        common = (general.get("mi_me"), general.get("dt_omega_pe"))
        labels = tuple(sorted(sections))
        if expected_general is None:
            expected_general = common
            expected_labels = labels
        elif common != expected_general or labels != expected_labels:
            raise SystemExit("Cases do not share mi/me, dt, and time sections")

        for label in labels:
            section = sections[label]
            records.append(
                {
                    "seed": seed,
                    "time_label": label,
                    "omega_pi_t": require_float(section, "omega_pi_t", summary_path),
                    "gamma_e": require_float(section, "gamma_e", summary_path),
                    "fit_standard_error": require_float(
                        section, "gamma_e_standard_error", summary_path
                    ),
                    "fit_r2": require_float(section, "fit_r2", summary_path),
                    "residual_standard_error": require_float(
                        section, "residual_standard_error", summary_path
                    ),
                    "fit_points": int(require_float(section, "fit_points", summary_path)),
                    "case_dir": str(case_dir),
                }
            )

    records.sort(key=lambda row: (row["omega_pi_t"], row["seed"]))
    return records, expected_general


def summarize(records):
    grouped = {}
    for row in records:
        grouped.setdefault(row["time_label"], []).append(row)

    summaries = []
    for label, rows in grouped.items():
        rows.sort(key=lambda row: row["seed"])
        values = [row["gamma_e"] for row in rows]
        if len(values) < 2:
            raise SystemExit("At least two seeds are required for a sample standard deviation")
        sample_sd = statistics.stdev(values)
        summaries.append(
            {
                "time_label": label,
                "omega_pi_t": rows[0]["omega_pi_t"],
                "n_seeds": len(rows),
                "gamma_mean": statistics.mean(values),
                "gamma_sample_standard_deviation": sample_sd,
                "gamma_standard_error_of_mean": sample_sd / math.sqrt(len(values)),
                "gamma_min": min(values),
                "gamma_max": max(values),
                "mean_fit_standard_error": statistics.mean(
                    row["fit_standard_error"] for row in rows
                ),
                "mean_fit_r2": statistics.mean(row["fit_r2"] for row in rows),
                "mean_residual_standard_error": statistics.mean(
                    row["residual_standard_error"] for row in rows
                ),
            }
        )
    summaries.sort(key=lambda row: row["omega_pi_t"])
    return summaries


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text_summary(path, records, summaries, general, args):
    seeds = sorted({row["seed"] for row in records})
    lines = [
        f"mi_me = {general[0]}",
        f"dt_omega_pe = {general[1]}",
        "seeds = " + ",".join(str(seed) for seed in seeds),
        f"n_seeds = {len(seeds)}",
        "gamma_error_definition = sample standard deviation across equal-weight seeds",
        "gamma_sem_definition = sample standard deviation / sqrt(n_seeds)",
    ]
    if args.reference is not None:
        lines.append(f"reference_gamma_e = {args.reference:.8g}")
    if args.reference_error is not None:
        lines.append(f"reference_gamma_e_error = {args.reference_error:.8g}")

    for summary in summaries:
        lines.extend(
            [
                "",
                f"[{summary['time_label']}]",
                f"omega_pi_t = {summary['omega_pi_t']:.8g}",
                f"gamma_e_mean = {summary['gamma_mean']:.8g}",
                "gamma_e_sample_standard_deviation = "
                f"{summary['gamma_sample_standard_deviation']:.8g}",
                "gamma_e_standard_error_of_mean = "
                f"{summary['gamma_standard_error_of_mean']:.8g}",
                f"gamma_e_min = {summary['gamma_min']:.8g}",
                f"gamma_e_max = {summary['gamma_max']:.8g}",
                f"mean_fit_standard_error = {summary['mean_fit_standard_error']:.8g}",
                f"mean_fit_r2 = {summary['mean_fit_r2']:.8g}",
                "mean_residual_standard_error = "
                f"{summary['mean_residual_standard_error']:.8g}",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_results(path, records, summaries, args):
    labels = [row["time_label"] for row in summaries]
    x_values = list(range(len(labels)))
    seeds = sorted({row["seed"] for row in records})

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for seed in seeds:
        seed_rows = {row["time_label"]: row for row in records if row["seed"] == seed}
        y_values = [seed_rows[label]["gamma_e"] for label in labels]
        ax.plot(x_values, y_values, marker="o", linewidth=1.2, alpha=0.75, label=f"seed {seed}")

    means = [row["gamma_mean"] for row in summaries]
    errors = [row["gamma_sample_standard_deviation"] for row in summaries]
    ax.errorbar(
        x_values,
        means,
        yerr=errors,
        color="black",
        marker="s",
        markersize=6,
        linewidth=2,
        capsize=5,
        label="mean +/- sample SD",
        zorder=5,
    )

    if args.reference is not None:
        ax.axhline(args.reference, color="tab:red", linestyle="--", linewidth=1.2, label="reference")
        if args.reference_error is not None:
            ax.axhspan(
                args.reference - args.reference_error,
                args.reference + args.reference_error,
                color="tab:red",
                alpha=0.10,
                linewidth=0,
            )

    ax.set_xticks(x_values, labels)
    ax.set_xlabel("Normalized time")
    ax.set_ylabel(r"Effective electron polytropic index $\gamma_e$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main():
    args = parse_args()
    if args.reference_error is not None and args.reference is None:
        raise SystemExit("--reference-error requires --reference")
    if args.reference_error is not None and args.reference_error < 0:
        raise SystemExit("--reference-error must be non-negative")

    records, general = load_cases(args.case_dirs)
    summaries = summarize(records)
    output_dir = args.output_dir or default_output_dir(args.case_dirs)
    output_dir.mkdir(parents=True, exist_ok=True)

    result_fields = [
        "seed",
        "time_label",
        "omega_pi_t",
        "gamma_e",
        "fit_standard_error",
        "fit_r2",
        "residual_standard_error",
        "fit_points",
        "case_dir",
    ]
    statistic_fields = list(summaries[0])
    write_csv(output_dir / "seed_results.csv", result_fields, records)
    write_csv(output_dir / "seed_statistics.csv", statistic_fields, summaries)
    write_text_summary(
        output_dir / "seed_statistics.txt", records, summaries, general, args
    )
    plot_results(output_dir / "seed_gamma.png", records, summaries, args)

    print(f"wrote {output_dir / 'seed_results.csv'}")
    print(f"wrote {output_dir / 'seed_statistics.csv'}")
    print(f"wrote {output_dir / 'seed_statistics.txt'}")
    print(f"wrote {output_dir / 'seed_gamma.png'}")


if __name__ == "__main__":
    main()
