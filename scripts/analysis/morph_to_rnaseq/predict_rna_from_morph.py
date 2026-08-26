#!/usr/bin/env python3
"""Fit/save morph→RNA bundle and predict from a morphology matrix only (no RNA)."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
from sklearn.pipeline import Pipeline

from run_morph_to_rnaseq import (
    TIME_ORDER,
    fit_cell_helper,
    image_rows,
    make_pipeline,
    parse_top_k_names,
    pooled_with_helper,
    read_cells,
    write_csv,
)


def load_stable_features(path: Path, min_folds: int = 4) -> list[str]:
    names = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(float(row["n_folds"])) >= min_folds:
                names.append(row["feature"])
    if not names:
        raise SystemExit(f"No features with n_folds>={min_folds} in {path}")
    return names


def fill_nan(x: np.ndarray) -> np.ndarray:
    out = np.array(x, dtype=float, copy=True)
    for j in range(out.shape[1]):
        med = np.nanmedian(out[:, j])
        out[~np.isfinite(out[:, j]), j] = med if np.isfinite(med) else 0.0
    return out


def pool_cells(
    morph_csv: Path,
    wanted_morph: list[str],
    cell_feats: list[str],
    drop_unknown_time: bool = True,
):
    features, by_feat, cell_images, cell_times, cell_clusters = read_cells(
        morph_csv, all_features=True, feature_names=cell_feats or None
    )
    names, x_img, y_img, image_ids, cluster_levels = image_rows(
        features,
        by_feat,
        cell_images,
        cell_times,
        cell_clusters,
        collinear_thresh=None,
        wanted_cols=wanted_morph,
        drop_unknown_time=drop_unknown_time,
    )
    cell_ok = np.array([t in TIME_ORDER for t in cell_times])
    cell_x = fill_nan(np.column_stack([by_feat[n] for n in features]))
    return {
        "features": features,
        "names": names,
        "x_img": x_img,
        "y_img": y_img,
        "image_ids": image_ids,
        "cluster_levels": cluster_levels,
        "cell_images": cell_images,
        "cell_times": cell_times,
        "cell_x": cell_x,
        "cell_ok": cell_ok,
    }


def assemble_x(pooled: dict, helper: Pipeline, helper_times: tuple[str, ...], model_cols: list[str]) -> np.ndarray:
    x = pooled["x_img"]
    names = list(pooled["names"])
    if helper is not None and helper_times:
        x = pooled_with_helper(
            x,
            pooled["image_ids"],
            pooled["cell_images"],
            pooled["cell_x"],
            helper,
            TIME_ORDER,
            helper_times,
        )
        names = names + [f"mean_Pcell_{t}" for t in helper_times]
    col_index = {n: i for i, n in enumerate(names)}
    missing = [n for n in model_cols if n not in col_index]
    if missing:
        raise SystemExit("Missing pooled columns: " + ", ".join(missing[:12]))
    return x[:, [col_index[n] for n in model_cols]]


def cmd_fit(args: argparse.Namespace) -> None:
    root: Path = args.project_root
    nested = args.nested_dir or (root / "results_nested_top_features")
    stable = load_stable_features(nested / "nested_feature_frequency.csv", args.min_folds)
    cell_feats, wanted_morph, helper_times = parse_top_k_names(stable)
    pooled = pool_cells(args.morph_csv, wanted_morph, cell_feats)
    y_idx = np.array([TIME_ORDER.index(t) for t in pooled["y_img"]])

    helper = None
    if helper_times:
        helper = fit_cell_helper(
            pooled["cell_x"][pooled["cell_ok"]],
            pooled["cell_times"][pooled["cell_ok"]],
            TIME_ORDER,
            C=1.0,
        )
    x = assemble_x(pooled, helper, helper_times, stable)
    clf = make_pipeline(args.C)
    clf.fit(x, y_idx)

    atlas = np.load(nested / "rna_by_time.npz", allow_pickle=True)
    bundle_dir = args.bundle_dir or (root / "bundle")
    bundle_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "clf": clf,
            "helper": helper,
            "helper_times": helper_times,
            "model_cols": stable,
            "wanted_morph": wanted_morph,
            "cell_feats": cell_feats,
            "times": list(TIME_ORDER),
            "note": "Input is morph only. RNA atlas is frozen from the BMP50 time course.",
        },
        bundle_dir / "model_bundle.joblib",
    )
    np.savez(
        bundle_dir / "rna_by_time.npz",
        times=atlas["times"],
        program_names=atlas["program_names"],
        programs=atlas["programs"],
        genes=atlas["genes"],
        gene_vst=atlas["gene_vst"],
    )
    (bundle_dir / "README.txt").write_text(
        "Predict RNA programs from morphology only.\n"
        "python predict_rna_from_morph.py predict --morph-csv CELLS.csv --bundle-dir this_dir\n"
        "No RNA-seq file is required.\n",
        encoding="utf-8",
    )
    print(json.dumps({"n_slices": int(len(y_idx)), "n_features": len(stable), "bundle": str(bundle_dir)}))


def cmd_predict(args: argparse.Namespace) -> None:
    bundle_dir: Path = args.bundle_dir
    bundle = joblib.load(bundle_dir / "model_bundle.joblib")
    atlas = np.load(bundle_dir / "rna_by_time.npz", allow_pickle=True)
    pooled = pool_cells(
        args.morph_csv,
        bundle["wanted_morph"],
        bundle["cell_feats"],
        drop_unknown_time=False,
    )
    x = assemble_x(pooled, bundle["helper"], tuple(bundle["helper_times"]), bundle["model_cols"])
    proba = bundle["clf"].predict_proba(x)
    pred_idx = np.argmax(proba, axis=1)
    programs = np.asarray(atlas["programs"], dtype=float)
    genes = [str(g) for g in atlas["genes"].tolist()]
    gene_vst = np.asarray(atlas["gene_vst"], dtype=float)
    program_names = [str(p) for p in atlas["program_names"].tolist()]
    pred_prog = proba @ programs
    pred_genes = proba @ gene_vst
    times = list(TIME_ORDER)
    true_times = [str(t) for t in pooled["y_img"]]
    sample_ids = [str(s) for s in pooled["image_ids"]]

    out_dir = args.out_dir or (bundle_dir / "predictions")
    write_csv(
        out_dir / "predicted_time.csv",
        ["slice", "sample", "true_time", "pred_time", "correct"]
        + [f"P_{t}" for t in times],
        [
            {
                "slice": i + 1,
                "sample": sample_ids[i],
                "true_time": true_times[i],
                "pred_time": times[pred_idx[i]],
                "correct": (
                    ""
                    if true_times[i] not in times
                    else str(true_times[i] == times[pred_idx[i]]).lower()
                ),
                **{f"P_{times[k]}": proba[i, k] for k in range(len(times))},
            }
            for i in range(len(pred_idx))
        ],
    )
    write_csv(
        out_dir / "predicted_programs.csv",
        ["slice", "sample", "pred_time"] + program_names,
        [
            {
                "slice": i + 1,
                "sample": sample_ids[i],
                "pred_time": times[pred_idx[i]],
                **{program_names[j]: pred_prog[i, j] for j in range(len(program_names))},
            }
            for i in range(len(pred_idx))
        ],
    )
    write_csv(
        out_dir / "predicted_panel_genes.csv",
        ["slice", "sample", "pred_time"] + genes,
        [
            {
                "slice": i + 1,
                "sample": sample_ids[i],
                "pred_time": times[pred_idx[i]],
                **{genes[j]: pred_genes[i, j] for j in range(len(genes))},
            }
            for i in range(len(pred_idx))
        ],
    )
    labeled = [i for i, t in enumerate(true_times) if t in times]
    n_correct = sum(1 for i in labeled if true_times[i] == times[pred_idx[i]])
    accuracy = (n_correct / len(labeled)) if labeled else None
    summary = {
        "n_slices": int(len(pred_idx)),
        "n_labeled": int(len(labeled)),
        "n_correct": int(n_correct) if labeled else None,
        "accuracy_labeled": accuracy,
        "out": str(out_dir),
    }
    (out_dir / "prediction_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary))
    for i in range(min(5, len(pred_idx))):
        print(
            f"{sample_ids[i]}: true={true_times[i]} pred={times[pred_idx[i]]}  "
            + " ".join(f"{t}={proba[i,k]:.2f}" for k, t in enumerate(times))
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    fit = sub.add_parser("fit")
    fit.add_argument("--project-root", type=Path, default=Path("/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq"))
    fit.add_argument(
        "--morph-csv",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results/single_cell_full_features.csv"),
    )
    fit.add_argument("--nested-dir", type=Path, default=None)
    fit.add_argument("--bundle-dir", type=Path, default=None)
    fit.add_argument("--min-folds", type=int, default=4)
    fit.add_argument("--C", type=float, default=0.1)
    fit.set_defaults(func=cmd_fit)

    pred = sub.add_parser("predict")
    pred.add_argument("--morph-csv", type=Path, required=True)
    pred.add_argument(
        "--bundle-dir",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/bundle"),
    )
    pred.add_argument("--out-dir", type=Path, default=None)
    pred.set_defaults(func=cmd_predict)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
