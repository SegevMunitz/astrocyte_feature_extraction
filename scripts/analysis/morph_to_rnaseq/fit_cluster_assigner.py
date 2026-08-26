#!/usr/bin/env python3
"""Fit a frozen cell→cluster assigner from an existing labeled morph CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SKIP = {
    "sample",
    "original_cell_id",
    "compact_cell_id",
    "Time",
    "cluster",
    "measurement_backend",
    "crop_path",
    "mask_path",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--morph-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--C", type=float, default=1.0)
    args = parser.parse_args()

    with args.morph_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit("Empty morph CSV")

    feature_names = [
        k
        for k in rows[0].keys()
        if k not in SKIP
        and (k.startswith("AreaShape_") or k.startswith("Intensity_") or k in {
            "skeleton_length_pixels",
            "nontrunk_branch_points",
        })
    ]
    x_rows = []
    y = []
    for row in rows:
        label = str(row.get("cluster", "")).strip()
        if not label:
            continue
        vals = []
        for name in feature_names:
            raw = row.get(name, "")
            try:
                vals.append(float(raw) if raw != "" else float("nan"))
            except ValueError:
                vals.append(float("nan"))
        x_rows.append(vals)
        y.append(label)

    if len(set(y)) < 2:
        raise SystemExit(f"Need ≥2 cluster labels; got {sorted(set(y))}")

    x = np.asarray(x_rows, dtype=float)
    for j in range(x.shape[1]):
        med = np.nanmedian(x[:, j])
        x[~np.isfinite(x[:, j]), j] = med if np.isfinite(med) else 0.0

    clf = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "logit",
                LogisticRegression(
                    max_iter=2000,
                    C=args.C,
                    solver="lbfgs",
                ),
            ),
        ]
    )
    clf.fit(x, y)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "clf": clf,
            "feature_names": feature_names,
            "classes": sorted(set(y), key=str),
            "note": "Assign morphology cluster labels for frac_cluster_* features.",
        },
        args.out,
    )
    print(
        json.dumps(
            {
                "n_cells": len(y),
                "n_features": len(feature_names),
                "classes": sorted(set(y), key=str),
                "out": str(args.out),
            }
        )
    )


if __name__ == "__main__":
    main()
