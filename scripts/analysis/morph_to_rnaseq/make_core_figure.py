#!/usr/bin/env python3
"""Build 6-panel core + supplementary RNA program figures for lab meeting."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_morph_to_rnaseq import TIME_ORDER, _f, plot_programs

CORE_PROGRAMS = (
    "BMP_Id",
    "cell_cycle",
    "reactive_GFAP",
    "cytoskeleton",
    "OPC_pdgfra",
    "late_7d",
)
SUPP_PROGRAMS = (
    "MHC_class_I",
    "Stat3_reactive",
    "glutamate_ion",
    "ECM_remodeling",
)


def read_program_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    names = [k for k in rows[0] if k != "time"]
    mat = np.array([[_f(row[n]) for n in names] for row in rows], dtype=float)
    return names, mat


def pred_by_time(nested: Path, program_names: list[str]) -> np.ndarray:
    cv_rows = list(csv.DictReader((nested / "cv_predictions.csv").open(encoding="utf-8")))
    pred_rows = list(csv.DictReader((nested / "predicted_programs.csv").open(encoding="utf-8")))
    pred_by_slice = {row["slice"]: row for row in pred_rows}
    out = []
    for t in TIME_ORDER:
        vals = [
            [_f(pred_by_slice[r["slice"]][p]) for p in program_names]
            for r in cv_rows
            if r["true_time"] == t
        ]
        out.append(np.nanmean(vals, axis=0) if vals else np.full(len(program_names), np.nan))
    return np.asarray(out, dtype=float)


def trajectory_notes(scores: np.ndarray, names: list[str]) -> list[dict[str, str]]:
    rows = []
    for j, name in enumerate(names):
        vals = scores[:, j]
        peak_i = int(np.nanargmax(vals))
        trough_i = int(np.nanargmin(vals))
        rows.append(
            {
                "program": name,
                "peak_time": TIME_ORDER[peak_i],
                "trough_time": TIME_ORDER[trough_i],
                "ctrl": f"{vals[0]:.2f}",
                "4h": f"{vals[1]:.2f}",
                "24h": f"{vals[2]:.2f}",
                "72h": f"{vals[3]:.2f}",
                "7d": f"{vals[4]:.2f}",
            }
        )
    return rows


def bar_correlations(metrics_path: Path, programs: tuple[str, ...], out: Path, title: str) -> None:
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    vals = [metrics["programs"][p]["pearson"] for p in programs]
    fig, ax = plt.subplots(figsize=(max(6.0, 0.65 * len(programs)), 3.5))
    ax.bar(list(programs), vals, color="#6a3d9a")
    ax.set_ylim(0, 1.05)
    ax.axhline(0.9, color="#cccccc", ls="--", lw=0.8)
    ax.set_ylabel("Pearson r")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--nested-dir",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/results_nested_top_features"),
    )
    args = parser.parse_args()
    nested = args.nested_dir
    plots = nested / "plots"
    names, actual = read_program_matrix(nested / "rna_program_scores.csv")
    name_index = {n: i for i, n in enumerate(names)}

    core_idx = [name_index[p] for p in CORE_PROGRAMS]
    actual_core = actual[:, core_idx]
    pred_core = pred_by_time(nested, list(CORE_PROGRAMS))
    plot_programs(
        actual_core,
        pred_core,
        list(CORE_PROGRAMS),
        plots / "core_six_programs_actual_vs_pred.png",
    )

    supp_idx = [name_index[p] for p in SUPP_PROGRAMS if p in name_index]
    supp_names = [names[i] for i in supp_idx]
    plot_programs(
        actual[:, supp_idx],
        pred_by_time(nested, supp_names),
        supp_names,
        plots / "supplementary_four_programs_actual_vs_pred.png",
    )

    bar_correlations(
        nested / "rna_metrics.json",
        CORE_PROGRAMS,
        plots / "core_six_correlations.png",
        "Core program recovery (morph → RNA)",
    )

    notes = trajectory_notes(actual, names)
    (nested / "program_trajectories.json").write_text(json.dumps(notes, indent=2), encoding="utf-8")
    print(json.dumps({"core_figure": str(plots / "core_six_programs_actual_vs_pred.png"), "trajectories": notes}, indent=2))


if __name__ == "__main__":
    main()
