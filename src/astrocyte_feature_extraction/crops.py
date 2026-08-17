"""Export isolated single-cell crops while preserving full-resolution masks."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from skimage.transform import resize

from .images import load_hwc
from .inventory import ImageRecord


class CropError(RuntimeError):
    """Raised when image and label geometry are inconsistent."""


def _write_channels(path: Path, array_hwc: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        np.moveaxis(array_hwc, -1, 0),
        metadata={"axes": "CYX"},
        photometric="minisblack",
    )


def _to_unit_interval(image: np.ndarray) -> np.ndarray:
    if np.issubdtype(image.dtype, np.integer):
        maximum = float(np.iinfo(image.dtype).max)
        return image.astype(np.float32) / maximum
    converted = image.astype(np.float32)
    if np.nanmin(converted) < 0 or np.nanmax(converted) > 1:
        raise CropError(
            "Article-compatible conversion requires floating-point images already scaled to [0, 1]"
        )
    return converted


def article_preprocess(
    crop_hwc: np.ndarray,
    size: int = 64,
    epsilon: float = 1.0e-5,
) -> np.ndarray:
    """Match the cited cell-image resize and per-channel z-normalization."""
    unit = _to_unit_interval(crop_hwc)
    quantized = np.rint(np.clip(unit, 0, 1) * 255).astype(np.uint8).astype(np.float32) / 255.0
    resized = resize(
        quantized,
        (size, size, quantized.shape[-1]),
        order=0,
        mode="constant",
        cval=0,
        anti_aliasing=False,
        preserve_range=True,
    ).astype(np.float32)
    means = resized.mean(axis=(0, 1), keepdims=True)
    standard_deviations = resized.std(axis=(0, 1), keepdims=True)
    return (resized - means) / (standard_deviations + float(epsilon))


def export_cell_crops(
    record: ImageRecord,
    labels_path: Path,
    output_root: Path,
    crop_config: dict[str, Any],
) -> pd.DataFrame:
    """Write cell-only native crops and return their global geometry."""
    image = load_hwc(record)
    labels = np.asarray(tifffile.imread(labels_path))
    if labels.shape != image.shape[:2]:
        raise CropError(
            f"Image/label shape mismatch for {record.image_id}: {image.shape[:2]} vs {labels.shape}"
        )
    margin = int(crop_config["margin_pixels"])
    if margin < 0:
        raise CropError("crops.margin_pixels must be non-negative")
    records: list[dict[str, Any]] = []
    for cell_id, location in enumerate(ndimage.find_objects(labels), start=1):
        if location is None:
            continue
        row_slice, column_slice = location
        y0 = max(0, int(row_slice.start) - margin)
        y1 = min(labels.shape[0], int(row_slice.stop) + margin)
        x0 = max(0, int(column_slice.start) - margin)
        x1 = min(labels.shape[1], int(column_slice.stop) + margin)
        local_labels = labels[y0:y1, x0:x1]
        selected = local_labels == cell_id
        if not selected.any():
            raise CropError(f"Cell {cell_id} disappeared during crop extraction")
        crop = image[y0:y1, x0:x1].copy()
        crop[~selected] = 0
        crop_dir = output_root / "cells" / record.image_id
        crop_path = crop_dir / f"cell_{cell_id:05d}.tif"
        mask_path = crop_dir / f"cell_{cell_id:05d}_mask.tif"
        _write_channels(crop_path, crop)
        tifffile.imwrite(
            mask_path,
            selected.astype(np.uint8),
            metadata={"axes": "YX"},
            photometric="minisblack",
        )
        article_path: str | None = None
        if crop_config.get("write_article_64x64", True):
            normalized = article_preprocess(
                crop,
                size=int(crop_config.get("article_size", 64)),
                epsilon=float(crop_config.get("article_epsilon", 1.0e-5)),
            )
            article_file = crop_dir / f"cell_{cell_id:05d}_article64.tif"
            _write_channels(article_file, normalized)
            article_path = str(article_file)
        records.append(
            {
                "image_id": record.image_id,
                "cell_id": cell_id,
                "crop_path": str(crop_path),
                "crop_mask_path": str(mask_path),
                "article_crop_path": article_path,
                "bbox_y0": y0,
                "bbox_x0": x0,
                "bbox_y1": y1,
                "bbox_x1": x1,
                "cell_area_pixels": int(selected.sum()),
            }
        )
    return pd.DataFrame.from_records(records)


def export_all_crops(
    records: list[ImageRecord],
    mask_paths: dict[str, Path],
    output_root: Path,
    crop_config: dict[str, Any],
) -> pd.DataFrame:
    """Export isolated cells for every segmented image."""
    tables = [
        export_cell_crops(record, mask_paths[record.image_id], output_root, crop_config)
        for record in records
    ]
    table = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    table_path = output_root / "tables" / "crops.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(table_path, index=False)
    return table
