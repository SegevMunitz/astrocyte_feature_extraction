"""Flatten pipeline cell crops for crop-level CellProfiler + skeleton measure."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .inventory import ImageRecord, discover_images
from .config import output_path

TIME_TOKEN = re.compile(r"(?<![A-Za-z0-9])(ctrl|4h|24h|72h|7d)(?![A-Za-z0-9])", re.IGNORECASE)


class CropMeasureError(RuntimeError):
    """Raised when crop flattening or measurement fails."""


def parse_time_label(name: str) -> str:
    match = TIME_TOKEN.search(name)
    return match.group(1).lower() if match else "unknown"


def sample_name_for_record(record: ImageRecord) -> str:
    return Path(record.relative_path).stem


def flatten_pipeline_crops(
    records: list[ImageRecord],
    output_root: Path,
    flat_crops_dir: Path,
    flat_masks_dir: Path,
) -> int:
    """Symlink cells/<image_id>/cell_*.tif into {sample}__cell_{id}.tif layout."""
    if flat_crops_dir.exists():
        shutil.rmtree(flat_crops_dir)
    if flat_masks_dir.exists():
        shutil.rmtree(flat_masks_dir)
    flat_crops_dir.mkdir(parents=True, exist_ok=True)
    flat_masks_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    for record in records:
        crop_dir = output_root / "cells" / record.image_id
        if not crop_dir.is_dir():
            raise CropMeasureError(f"Missing crop directory: {crop_dir}")
        sample = sample_name_for_record(record)
        for crop_path in sorted(crop_dir.glob("cell_*.tif")):
            if crop_path.name.endswith("_article64.tif") or crop_path.name.endswith("_mask.tif"):
                continue
            stem = crop_path.stem  # cell_00001
            cell_id = int(stem.split("_")[-1])
            mask_path = crop_dir / f"{stem}_mask.tif"
            if not mask_path.is_file():
                raise CropMeasureError(f"Missing mask for {crop_path}")
            out_crop = flat_crops_dir / f"{sample}__cell_{cell_id:05d}.tif"
            out_mask = flat_masks_dir / f"{sample}__cell_{cell_id:05d}_mask.tif"
            out_crop.symlink_to(crop_path.resolve())
            out_mask.symlink_to(mask_path.resolve())
            n += 1
    return n


def apply_cluster_assigner(morph_csv: Path, assigner_path: Path) -> None:
    """Fill missing cluster labels using a frozen sklearn assigner."""
    import csv

    import joblib
    import numpy as np

    payload = joblib.load(assigner_path)
    feature_names: list[str] = list(payload["feature_names"])
    clf = payload["clf"]
    classes = [str(c) for c in payload["classes"]]

    with morph_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        if not rows:
            return
        fieldnames = list(rows[0].keys())
        if "cluster" not in fieldnames:
            fieldnames.insert(fieldnames.index("Time") + 1 if "Time" in fieldnames else 0, "cluster")

    matrix = []
    need_idx = []
    for i, row in enumerate(rows):
        if str(row.get("cluster", "")).strip():
            continue
        need_idx.append(i)
        vals = []
        for name in feature_names:
            raw = row.get(name, "")
            try:
                vals.append(float(raw) if raw != "" else float("nan"))
            except ValueError:
                vals.append(float("nan"))
        matrix.append(vals)

    if need_idx:
        x = np.asarray(matrix, dtype=float)
        for j in range(x.shape[1]):
            med = np.nanmedian(x[:, j])
            x[~np.isfinite(x[:, j]), j] = med if np.isfinite(med) else 0.0
        pred = clf.predict(x)
        for row_i, label in zip(need_idx, pred):
            rows[row_i]["cluster"] = str(label)

    # Normalize Time parse for BMP4_* names if blank/unknown
    for row in rows:
        sample = str(row.get("sample", ""))
        time = str(row.get("Time", "")).strip()
        if not time or time == "unknown":
            row["Time"] = parse_time_label(sample)

    with morph_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    # silence unused
    _ = classes


def run_crop_measure(config: dict[str, Any]) -> Path:
    """Flatten crops, run crop-level CP+skeleton measure, assign clusters."""
    records = discover_images(config)
    output_root = output_path(config)
    flat_crops = output_root / "tables" / "flat_crops"
    flat_masks = output_root / "tables" / "flat_masks"
    morph_csv = output_root / "tables" / "single_cell_full_features.csv"

    n = flatten_pipeline_crops(records, output_root, flat_crops, flat_masks)
    print(json.dumps({"flattened_cells": n, "crops": str(flat_crops)}))

    cp = config["cellprofiler"]
    measure_script = Path(config["paths"].get(
        "crop_measure_script",
        Path(__file__).resolve().parents[2] / "scripts" / "analysis" / "measure_single_cell_crops_full.py",
    ))
    python_exe = Path(cp.get("python_executable") or cp["executable"]).expanduser()
    # Prefer dedicated CP python; executable may be cellprofiler binary
    if python_exe.name == "cellprofiler":
        python_exe = python_exe.parent / "python"

    workers = int(config.get("crop_measure", {}).get("workers", 8))
    cmd = [
        str(python_exe),
        str(measure_script),
        "--crops-dir",
        str(flat_crops),
        "--masks-dir",
        str(flat_masks),
        "--output",
        str(morph_csv),
        "--workers",
        str(workers),
    ]
    existing = config["paths"].get("cluster_csv")
    if existing:
        cmd.extend(["--existing-features", str(Path(existing).expanduser())])

    print(json.dumps({"crop_measure_cmd": cmd}))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise CropMeasureError(f"crop measure failed with exit code {result.returncode}")

    assigner = config["paths"].get("cluster_assigner")
    if assigner and Path(assigner).expanduser().is_file():
        apply_cluster_assigner(morph_csv, Path(assigner).expanduser())
        print(json.dumps({"cluster_assigner": str(assigner)}))

    return morph_csv


def run_predict_rna(config: dict[str, Any]) -> Path:
    """Predict RNA programs for all FOVs from morph CSV + frozen bundle."""
    output_root = output_path(config)
    morph_csv = Path(
        config["paths"].get("morph_csv")
        or (output_root / "tables" / "single_cell_full_features.csv")
    ).expanduser()
    if not morph_csv.is_file():
        raise CropMeasureError(f"Missing morph CSV: {morph_csv}")
    bundle_dir = Path(config["paths"]["bundle_dir"]).expanduser()
    out_dir = output_root / "rna_predictions"
    out_dir.mkdir(parents=True, exist_ok=True)

    predict_script = Path(__file__).resolve().parents[2] / "scripts" / "analysis" / "morph_to_rnaseq" / "predict_rna_from_morph.py"
    # Use pipeline venv python (sklearn/joblib)
    import sys

    cmd = [
        sys.executable,
        str(predict_script),
        "predict",
        "--morph-csv",
        str(morph_csv),
        "--bundle-dir",
        str(bundle_dir),
        "--out-dir",
        str(out_dir),
    ]
    print(json.dumps({"predict_cmd": cmd}))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise CropMeasureError(f"predict_rna failed with exit code {result.returncode}")
    return out_dir
