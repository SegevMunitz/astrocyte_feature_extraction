#!/usr/bin/env python3
"""Score nested morph→RNA programs and export frozen gene-level atlas tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from run_morph_to_rnaseq import (
    TIME_ORDER,
    TIME_TO_RNA,
    assess_programs,
    infer_groups_from_names,
    load_program_defs,
    plot_programs,
    program_matrix,
    read_col_groups,
    read_vst,
    resolve_programs_path,
    write_csv,
    _f,
)

CONTRAST_FILES = {
    "4h": "4H_vs_ctrl.csv",
    "24h": "24H_vs_ctrl.csv",
    "72h": "72H_vs_ctrl.csv",
    "7d": "7D_vs_ctrl.csv",
}


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 5:
        return float("nan")
    a = a[mask]
    b = b[mask]
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def resolve_groups(samples: list[str], coldata_csv: Path) -> dict[str, str]:
    groups = read_col_groups(coldata_csv, samples)
    inferred = infer_groups_from_names(samples)
    for sample in samples:
        if sample not in groups or groups[sample] not in TIME_TO_RNA.values():
            if inferred.get(sample):
                groups[sample] = inferred[sample]
            elif inferred.get(sample.lstrip("X")):
                groups[sample] = inferred[sample.lstrip("X")]
    return groups


def panel_genes(programs: dict[str, tuple[str, ...]]) -> list[str]:
    seen: list[str] = []
    for genes in programs.values():
        for gene in genes:
            if gene not in seen:
                seen.append(gene)
    return seen


def load_contrast_genes(contrast_dir: Path, n_per_time: int = 8) -> list[str]:
    extra: list[str] = []
    if not contrast_dir.is_dir():
        return extra
    for time, fname in CONTRAST_FILES.items():
        path = contrast_dir / fname
        if not path.is_file():
            continue
        scored = []
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                gene = row.get("gene_id") or row.get("gene") or ""
                try:
                    padj = float(row["padj"])
                    lfc = float(row["log2FoldChange"])
                except (KeyError, TypeError, ValueError):
                    continue
                if math.isnan(padj) or math.isnan(lfc) or padj >= 0.05 or abs(lfc) < 1:
                    continue
                scored.append((abs(lfc), gene))
        scored.sort(reverse=True)
        for _abs, gene in scored[:n_per_time]:
            if gene and gene not in extra:
                extra.append(gene)
    return extra


def mean_vst_by_time(
    genes: list[str],
    samples: list[str],
    vst: np.ndarray,
    groups: dict[str, str],
    wanted: list[str],
) -> tuple[list[str], np.ndarray]:
    gene_index = {g: i for i, g in enumerate(genes)}
    keep = [g for g in wanted if g in gene_index]
    mat = np.zeros((len(TIME_ORDER), len(keep)), dtype=float)
    for t_i, time in enumerate(TIME_ORDER):
        rna_group = TIME_TO_RNA[time]
        cols = [
            j
            for j, sample in enumerate(samples)
            if groups.get(sample) == rna_group or groups.get(sample.lstrip("X")) == rna_group
        ]
        if not cols:
            mat[t_i, :] = np.nan
            continue
        for g_i, gene in enumerate(keep):
            mat[t_i, g_i] = float(np.nanmean(vst[gene_index[gene], cols]))
    return keep, mat


def plot_program_r(metrics: dict[str, dict[str, float]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(metrics)
    vals = [metrics[n]["pearson"] for n in names]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.7 * len(names)), 3.8))
    ax.bar(names, vals, color="#6a3d9a")
    ax.set_ylim(-1.05, 1.05)
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.set_ylabel("Pearson r (pred vs true-time program)")
    ax.set_title("RNA program recovery from morphology")
    ax.tick_params(axis="x", rotation=40, labelsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq"),
    )
    parser.add_argument("--nested-dir", type=Path, default=None)
    parser.add_argument("--vst-csv", type=Path, default=None)
    parser.add_argument("--coldata-csv", type=Path, default=None)
    parser.add_argument(
        "--contrast-dir",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2/results/contrasts"),
    )
    parser.add_argument("--programs-yaml", type=Path, default=None)
    args = parser.parse_args()
    root = args.project_root
    nested = args.nested_dir or (root / "results_nested_top_features")
    vst_csv = args.vst_csv or (root / "data" / "vst_matrix.csv")
    coldata_csv = args.coldata_csv or (root / "data" / "colData.csv")

    cv_rows = read_csv_rows(nested / "cv_predictions.csv")
    p_cols = [f"P_{t}" for t in TIME_ORDER]
    slices = [row["slice"] for row in cv_rows]
    true_times = [row["true_time"] for row in cv_rows]
    proba = np.array([[_f(row[c]) for c in p_cols] for row in cv_rows], dtype=float)

    genes, samples, vst = read_vst(vst_csv)
    groups = resolve_groups(samples, coldata_csv)
    program_defs = load_program_defs(args.programs_yaml, project_root=root)
    shutil.copy2(
        resolve_programs_path(args.programs_yaml, project_root=root),
        nested / "programs.yaml",
    )
    selected_programs, viability = assess_programs(
        program_defs, genes, samples, vst, groups
    )
    program_names, y_program, used_genes = program_matrix(
        genes, samples, vst, groups, programs=selected_programs
    )
    time_index = {t: i for i, t in enumerate(TIME_ORDER)}
    pred = proba @ y_program
    true = np.array([y_program[time_index[t]] for t in true_times], dtype=float)

    write_csv(
        nested / "rna_program_scores.csv",
        ["time"] + program_names,
        [
            {
                "time": TIME_ORDER[i],
                **{program_names[j]: y_program[i, j] for j in range(len(program_names))},
            }
            for i in range(len(TIME_ORDER))
        ],
    )
    write_csv(
        nested / "predicted_programs.csv",
        ["slice", "true_time"] + program_names,
        [
            {
                "slice": s,
                "true_time": t,
                **{program_names[j]: pred[i, j] for j in range(len(program_names))},
            }
            for i, (s, t) in enumerate(zip(slices, true_times))
        ],
    )

    per_program = {}
    for j, name in enumerate(program_names):
        r = _pearson(pred[:, j], true[:, j])
        mae = float(np.nanmean(np.abs(pred[:, j] - true[:, j])))
        per_program[name] = {"pearson": r, "mae": mae}
    dropped = [row["program"] for row in viability if not row["kept"]]
    rna_metrics = {
        "n_slices": len(cv_rows),
        "n_programs": len(program_names),
        "programs": per_program,
        "mean_pearson": float(np.nanmean([p["pearson"] for p in per_program.values()])),
        "program_genes_used": used_genes,
        "programs_dropped": dropped,
        "note": "Predicted program = P(time|morph) mixed with frozen bulk VST program scores. Morph→time was not retrained.",
    }
    (nested / "rna_metrics.json").write_text(json.dumps(rna_metrics, indent=2), encoding="utf-8")
    (nested / "program_viability.json").write_text(
        json.dumps(viability, indent=2, default=str), encoding="utf-8"
    )
    plot_program_r(per_program, nested / "plots" / "rna_program_correlations.png")
    pred_by_time = np.vstack(
        [
            pred[np.array(true_times) == t].mean(axis=0)
            if any(tt == t for tt in true_times)
            else np.full(len(program_names), np.nan)
            for t in TIME_ORDER
        ]
    )
    plot_programs(
        y_program, pred_by_time, program_names, nested / "plots" / "programs_actual_vs_pred.png"
    )

    wanted = panel_genes(selected_programs) + load_contrast_genes(args.contrast_dir)
    uniq: list[str] = []
    for gene in wanted:
        if gene not in uniq:
            uniq.append(gene)
    keep, gene_by_time = mean_vst_by_time(genes, samples, vst, groups, uniq)

    write_csv(
        nested / "rna_gene_atlas.csv",
        ["time"] + keep,
        [
            {"time": TIME_ORDER[i], **{keep[j]: gene_by_time[i, j] for j in range(len(keep))}}
            for i in range(len(TIME_ORDER))
        ],
    )

    pred_genes = proba @ gene_by_time
    ctrl = gene_by_time[0]
    write_csv(
        nested / "predicted_panel_genes.csv",
        ["slice", "true_time"] + keep,
        [
            {
                "slice": s,
                "true_time": t,
                **{keep[j]: pred_genes[i, j] for j in range(len(keep))},
            }
            for i, (s, t) in enumerate(zip(slices, true_times))
        ],
    )

    example_rows = []
    seen_times = set()
    for i, (s, t) in enumerate(zip(slices, true_times)):
        if t in seen_times:
            continue
        seen_times.add(t)
        delta = pred_genes[i] - ctrl
        order = np.argsort(np.abs(delta))[::-1][:12]
        for gene_i in order:
            example_rows.append(
                {
                    "slice": s,
                    "true_time": t,
                    "gene": keep[gene_i],
                    "predicted_vst": pred_genes[i, gene_i],
                    "delta_vs_ctrl": delta[gene_i],
                }
            )
    write_csv(
        nested / "example_fov_top_genes.csv",
        ["slice", "true_time", "gene", "predicted_vst", "delta_vs_ctrl"],
        example_rows,
    )

    np.savez(
        nested / "rna_by_time.npz",
        times=np.array(TIME_ORDER),
        program_names=np.array(program_names),
        programs=y_program,
        genes=np.array(keep),
        gene_vst=gene_by_time,
    )
    print(
        json.dumps(
            {
                "mean_pearson": rna_metrics["mean_pearson"],
                "n_programs": len(program_names),
                "dropped": dropped,
                "programs": per_program,
            },
            indent=2,
        )
    )
    print("Wrote RNA metrics and gene tables to", nested)


if __name__ == "__main__":
    main()
