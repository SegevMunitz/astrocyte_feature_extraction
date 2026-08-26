#!/usr/bin/env python3
"""True-time cytoskeleton RNA vs skeleton morph — time-course (non-linear)."""

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
from scipy import stats

from run_morph_to_rnaseq import TIME_ORDER, _f

SKELETON_FEATURES = (
    "median_skeleton_length_pixels",
    "p90_skeleton_length_pixels",
    "median_nontrunk_branch_points",
    "p90_nontrunk_branch_points",
)

TIME_COLORS = {
    "ctrl": "#4c4c4c",
    "4h": "#1f77b4",
    "24h": "#ff7f0e",
    "72h": "#2ca02c",
    "7d": "#d62728",
}


def stars(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def short_label(feature: str) -> str:
    return (
        feature.replace("median_", "med ")
        .replace("p90_", "p90 ")
        .replace("skeleton_length_pixels", "skel. length")
        .replace("nontrunk_branch_points", "branch pts")
    )


def mean_sem_by_time(values: np.ndarray, times: list[str]) -> tuple[np.ndarray, np.ndarray]:
    means = []
    sems = []
    for t in TIME_ORDER:
        mask = np.array([tt == t for tt in times])
        vals = values[mask]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            means.append(float("nan"))
            sems.append(float("nan"))
        else:
            means.append(float(np.mean(vals)))
            sems.append(float(stats.sem(vals)) if vals.size > 1 else 0.0)
    return np.asarray(means, dtype=float), np.asarray(sems, dtype=float)


def kruskal_by_time(values: np.ndarray, times: list[str]) -> tuple[float, float]:
    groups = []
    for t in TIME_ORDER:
        vals = values[np.array([tt == t for tt in times])]
        vals = vals[np.isfinite(vals)]
        if vals.size:
            groups.append(vals)
    if len(groups) < 2:
        return float("nan"), float("nan")
    h, p = stats.kruskal(*groups)
    return float(h), float(p)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--nested-dir",
        type=Path,
        default=Path(
            "/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/results_nested_top_features"
        ),
    )
    args = parser.parse_args()
    nested = args.nested_dir

    feat_rows = read_csv_rows(nested / "image_features.csv")
    pred_rows = read_csv_rows(nested / "predicted_programs.csv")
    score_rows = read_csv_rows(nested / "rna_program_scores.csv")
    true_by_time = {row["time"]: _f(row["cytoskeleton"]) for row in score_rows}
    time_by_slice = {row["slice"]: row["true_time"] for row in pred_rows}

    missing = [c for c in SKELETON_FEATURES if c not in feat_rows[0]]
    if missing:
        raise SystemExit(f"Missing skeleton columns in image_features.csv: {missing}")

    times: list[str] = []
    skeleton = {c: [] for c in SKELETON_FEATURES}
    for row in feat_rows:
        s = row["slice"]
        t = time_by_slice.get(s) or row.get("time", "")
        if t not in true_by_time:
            continue
        times.append(t)
        for c in SKELETON_FEATURES:
            skeleton[c].append(_f(row[c]))

    rna_by_time = np.array([true_by_time[t] for t in TIME_ORDER], dtype=float)
    rna_peak = TIME_ORDER[int(np.nanargmax(rna_by_time))]

    records = []
    x_idx = np.arange(len(TIME_ORDER))

    # Figure: 2x2 time courses — skeleton mean±SEM + true cytoskeleton on twin axis
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
    flat = axes.ravel()
    for i, feat in enumerate(SKELETON_FEATURES):
        ax = flat[i]
        vals = np.asarray(skeleton[feat], dtype=float)
        means, sems = mean_sem_by_time(vals, times)
        h, p = kruskal_by_time(vals, times)
        star = stars(p)
        skel_peak = TIME_ORDER[int(np.nanargmax(means))] if np.isfinite(means).any() else "?"

        ax.errorbar(
            x_idx,
            means,
            yerr=sems,
            fmt="-o",
            color="#333333",
            ecolor="#888888",
            capsize=3,
            lw=1.6,
            label="skeleton morph (mean±SEM)",
        )
        ax.set_ylabel(short_label(feat), fontsize=9, color="#333333")
        ax.tick_params(axis="y", labelcolor="#333333", labelsize=8)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(TIME_ORDER)
        ax.set_title(
            f"{short_label(feat)}  ·  peak {skel_peak}  ·  Kruskal {star}",
            fontsize=9,
        )

        ax2 = ax.twinx()
        ax2.plot(
            x_idx,
            rna_by_time,
            "--s",
            color="#6a3d9a",
            lw=1.8,
            markersize=6,
            label="true cytoskeleton RNA",
        )
        ax2.set_ylabel("true cytoskeleton (z)", fontsize=8, color="#6a3d9a")
        ax2.tick_params(axis="y", labelcolor="#6a3d9a", labelsize=8)
        ax2.axhline(0, color="#cccccc", lw=0.8)

        records.append(
            {
                "feature": feat,
                "rna_source": "true_time",
                "skel_means_by_time": {t: float(means[j]) for j, t in enumerate(TIME_ORDER)},
                "skel_sem_by_time": {t: float(sems[j]) for j, t in enumerate(TIME_ORDER)},
                "skel_peak_time": skel_peak,
                "rna_peak_time": rna_peak,
                "kruskal_H": h,
                "kruskal_p": p,
                "stars": star,
                "n_slices": int(np.isfinite(vals).sum()),
            }
        )

    # Shared legend
    handles = [
        plt.Line2D([0], [0], color="#333333", marker="o", lw=1.6, label="skeleton morph (mean±SEM)"),
        plt.Line2D(
            [0],
            [0],
            color="#6a3d9a",
            marker="s",
            ls="--",
            lw=1.8,
            label="true bulk cytoskeleton RNA (z)",
        ),
    ]
    fig.suptitle(
        "True cytoskeleton RNA vs skeleton morph over BMP50 time\n"
        "Black: skeleton mean±SEM · Purple: true bulk cytoskeleton (z) · "
        "Kruskal-Wallis across times: * p<0.05 ** p<0.01 *** p<0.001",
        fontsize=11,
        y=0.99,
    )
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))

    plots = nested / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    out_png = plots / "cytoskeleton_vs_skeleton_corr.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Second figure: scatter vs true RNA with quadratic (captures mid-peak), colored by time
    fig2, axes2 = plt.subplots(2, 2, figsize=(10.5, 7.0))
    flat2 = axes2.ravel()
    true_y = np.array([true_by_time[t] for t in times], dtype=float)
    colors = [TIME_COLORS.get(t, "#888888") for t in times]
    for i, feat in enumerate(SKELETON_FEATURES):
        ax = flat2[i]
        x = np.asarray(skeleton[feat], dtype=float)
        ax.scatter(x, true_y, c=colors, s=28, alpha=0.85, edgecolors="none")
        mask = np.isfinite(x) & np.isfinite(true_y)
        if mask.sum() >= 4 and np.nanstd(x[mask]) > 0:
            coef = np.polyfit(x[mask], true_y[mask], 2)
            xs = np.linspace(np.nanmin(x[mask]), np.nanmax(x[mask]), 80)
            ax.plot(xs, np.polyval(coef, xs), color="#6a3d9a", lw=1.6, label="quadratic fit")
            # R² of quadratic
            yhat = np.polyval(coef, x[mask])
            ss_res = float(np.sum((true_y[mask] - yhat) ** 2))
            ss_tot = float(np.sum((true_y[mask] - np.mean(true_y[mask])) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
            records[i]["quadratic_r2"] = r2
            ax.set_title(f"{short_label(feat)}\nquadratic R² = {r2:.2f}", fontsize=9)
        else:
            ax.set_title(short_label(feat), fontsize=9)
        ax.set_xlabel("skeleton morph", fontsize=8)
        if i % 2 == 0:
            ax.set_ylabel("true cytoskeleton (z)", fontsize=9)
        ax.tick_params(labelsize=8)

    time_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=TIME_COLORS[t],
            markersize=7,
            label=t,
        )
        for t in TIME_ORDER
    ]
    fig2.suptitle(
        "True cytoskeleton vs skeleton (per slice) — quadratic fit (not linear)\n"
        "Point color = true time; mid-course RNA peak can look like an inverted U",
        fontsize=11,
        y=0.99,
    )
    fig2.legend(
        handles=time_handles,
        loc="lower center",
        ncol=5,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
        title="time",
    )
    fig2.tight_layout(rect=(0, 0.06, 1, 0.90))
    out_scatter = plots / "cytoskeleton_vs_skeleton_scatter_quadratic.png"
    fig2.savefig(out_scatter, dpi=150, bbox_inches="tight")
    plt.close(fig2)

    out_json = nested / "cytoskeleton_skeleton_corr.json"
    payload = {
        "note": (
            "RNA is TRUE-TIME bulk cytoskeleton (not morph-predicted). "
            "Association is time-structured (often rise to ~24h then fall), so linear/Spearman "
            "slice correlations are secondary. Primary view = time courses + Kruskal-Wallis "
            "across times; scatter uses a quadratic fit for the mid-peak shape."
        ),
        "rna_source": "true_time",
        "rna_by_time": {t: float(true_by_time[t]) for t in TIME_ORDER},
        "rna_peak_time": rna_peak,
        "features": records,
    }
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "png_timecourse": str(out_png),
                "png_scatter": str(out_scatter),
                "json": str(out_json),
                "rna": "true_time",
                "rna_peak": rna_peak,
                "skel_peaks": {r["feature"]: r["skel_peak_time"] for r in records},
                "kruskal": {
                    r["feature"]: f"{r['stars']} (p={r['kruskal_p']:.3g})" for r in records
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
