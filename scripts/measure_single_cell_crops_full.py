"""Measure full CellProfiler SizeShape+Intensity features on per-cell crop TIFFs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from cellprofiler.modules.measureobjectintensity import MeasureObjectIntensity
from cellprofiler.modules.measureobjectsizeshape import MeasureObjectSizeShape
from cellprofiler_core.image import Image, ImageSetList
from cellprofiler_core.measurement import Measurements
from cellprofiler_core.object import ObjectSet, Objects
from cellprofiler_core.pipeline import Pipeline
from cellprofiler_core.preferences import set_headless
from cellprofiler_core.workspace import Workspace


def _load_image(path: Path) -> np.ndarray:
    array = np.asarray(tifffile.imread(path))
    if array.ndim == 2:
        return array[..., np.newaxis]
    if array.ndim == 3 and array.shape[0] <= 8 and array.shape[0] < array.shape[-1]:
        return np.moveaxis(array, 0, -1)
    if array.ndim == 3:
        return array
    raise ValueError(f"Unsupported crop shape {array.shape} at {path}")


def _scaled_plane(plane: np.ndarray) -> np.ndarray:
    if np.issubdtype(plane.dtype, np.integer):
        return plane.astype(np.float32) / float(np.iinfo(plane.dtype).max)
    result = plane.astype(np.float32)
    if np.nanmin(result) < 0 or np.nanmax(result) > 1:
        maximum = float(np.nanmax(result))
        if maximum <= 0:
            raise ValueError("Image plane has non-positive maximum")
        return result / maximum
    return result


def _labels_from_mask(mask: np.ndarray) -> np.ndarray:
    labels = np.asarray(mask)
    if labels.ndim != 2:
        raise ValueError(f"Expected 2D mask, got {labels.shape}")
    if labels.max() <= 1:
        return (labels > 0).astype(np.int32)
    return labels.astype(np.int32, copy=False)


def measure_crop(crop_path: Path, mask_path: Path) -> dict[str, float]:
    image_hwc = _load_image(crop_path)
    labels = _labels_from_mask(np.asarray(tifffile.imread(mask_path)))
    if labels.shape != image_hwc.shape[:2]:
        raise ValueError(
            f"Crop/mask shape mismatch for {crop_path.name}: "
            f"{image_hwc.shape[:2]} vs {labels.shape}"
        )

    pipeline = Pipeline()
    image_set_list = ImageSetList()
    image_set = image_set_list.get_image_set(0)
    channel_names: list[str] = []
    for channel_index in range(image_hwc.shape[-1]):
        name = f"Channel_{channel_index}"
        channel_names.append(name)
        image_set.add(
            name,
            Image(_scaled_plane(image_hwc[..., channel_index]), dimensions=2),
        )

    objects = Objects()
    objects.segmented = labels
    object_set = ObjectSet()
    object_set.add_objects(objects, "Cells")

    size_shape = MeasureObjectSizeShape()
    size_shape.set_module_num(1)
    size_shape.objects_list.value = "Cells"
    size_shape.calculate_zernikes.value = True
    size_shape.calculate_advanced.value = True
    pipeline.add_module(size_shape)

    intensity = MeasureObjectIntensity()
    intensity.set_module_num(2)
    intensity.images_list.value = channel_names
    intensity.objects_list.value = "Cells"
    pipeline.add_module(intensity)

    measurements = Measurements()
    for module in pipeline.modules():
        workspace = Workspace(
            pipeline,
            module,
            image_set,
            object_set,
            measurements,
            image_set_list,
        )
        module.run(workspace)

    row: dict[str, float] = {}
    for feature in measurements.get_feature_names("Cells"):
        values = np.asarray(measurements.get_measurement("Cells", feature))
        if len(values) != 1:
            raise ValueError(f"Expected one object in {crop_path.name}, got {len(values)}")
        row[feature] = float(values[0])
    return row


def _pair_crops_and_masks(crops_dir: Path, masks_dir: Path) -> list[tuple[Path, Path, str]]:
    pairs: list[tuple[Path, Path, str]] = []
    for crop_path in sorted(crops_dir.glob("*.tif*")):
        stem = crop_path.stem
        mask_candidates = [
            masks_dir / f"{stem}.tif",
            masks_dir / f"{stem}.tiff",
            masks_dir / f"{stem}_mask.tif",
            masks_dir / f"{stem}_mask.tiff",
        ]
        mask_path = next((path for path in mask_candidates if path.is_file()), None)
        if mask_path is None:
            raise FileNotFoundError(f"No mask found for crop {crop_path}")
        sample = stem
        pairs.append((crop_path, mask_path, sample))
    if not pairs:
        raise FileNotFoundError(f"No crop TIFFs found in {crops_dir}")
    return pairs


def _parse_time(sample: str) -> str | None:
    match = re.match(r"^(ctrl|4h|24h|72h|7d)", sample)
    return match.group(1) if match else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=Path(
            "/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis"
        ),
    )
    parser.add_argument(
        "--crops-dir",
        type=Path,
        default=None,
        help="Defaults to <analysis-root>/results/single_cell_crops",
    )
    parser.add_argument(
        "--masks-dir",
        type=Path,
        default=None,
        help="Defaults to <analysis-root>/results/single_cell_masks",
    )
    parser.add_argument(
        "--existing-features",
        type=Path,
        default=None,
        help="Optional prior CSV to merge cluster/metadata columns from",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Defaults to <analysis-root>/results/single_cell_full_features.csv",
    )
    parser.add_argument("--limit", type=int, default=None, help="Process first N cells")
    args = parser.parse_args()

    analysis_root = args.analysis_root.expanduser()
    crops_dir = args.crops_dir or (analysis_root / "results" / "single_cell_crops")
    masks_dir = args.masks_dir or (analysis_root / "results" / "single_cell_masks")
    output_path = args.output or (analysis_root / "results" / "single_cell_full_features.csv")

    set_headless()
    pairs = _pair_crops_and_masks(crops_dir, masks_dir)
    if args.limit is not None:
        pairs = pairs[: args.limit]

    existing: pd.DataFrame | None = None
    if args.existing_features is not None:
        existing_path = args.existing_features.expanduser()
        if existing_path.is_file():
            existing = pd.read_csv(existing_path)

    rows: list[dict[str, object]] = []
    feature_order: list[str] | None = None
    for index, (crop_path, mask_path, sample) in enumerate(pairs, start=1):
        measurements = measure_crop(crop_path, mask_path)
        if feature_order is None:
            feature_order = sorted(measurements)
        elif sorted(measurements) != feature_order:
            raise ValueError(f"Feature schema changed at {crop_path.name}")

        row: dict[str, object] = {
            "sample": sample,
            "crop_path": str(crop_path),
            "mask_path": str(mask_path),
            "Time": _parse_time(sample),
            **measurements,
        }
        if existing is not None and "sample" in existing.columns:
            match = existing.loc[existing["sample"] == sample]
            if len(match) == 1:
                for column in ("cluster", "original_cell_id", "compact_cell_id"):
                    if column in match.columns:
                        row[column] = match.iloc[0][column]
        rows.append(row)
        if index % 250 == 0:
            print(f"Measured {index}/{len(pairs)} cells")

    table = pd.DataFrame.from_records(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    feature_count = len(feature_order or [])
    print(f"Wrote {len(table)} rows x {feature_count} CellProfiler features to {output_path}")


if __name__ == "__main__":
    main()
