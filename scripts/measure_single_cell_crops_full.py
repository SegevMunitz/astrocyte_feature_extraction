"""Measure full CellProfiler SizeShape+Intensity features on per-cell crop TIFFs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


def _require_cellprofiler() -> None:
    try:
        import cellprofiler_core  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "CellProfiler Python packages are not importable. "
            "Activate the analysis venv that includes CellProfiler 4.2.8, e.g.\n"
            "  source .venv_astro_pipeline/bin/activate\n"
            "  python -c \"import cellprofiler; print(cellprofiler.__version__)\""
        ) from exc


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
    foreground = labels > 0
    if not foreground.any():
        raise ValueError("Mask has no foreground pixels")
    if labels.max() <= 1:
        return foreground.astype(np.int32)
    unique, counts = np.unique(labels[foreground], return_counts=True)
    largest = int(unique[int(np.argmax(counts))])
    return np.where(labels == largest, 1, 0).astype(np.int32)


def _stem_aliases(stem: str) -> list[str]:
    aliases = [stem]
    for suffix in ("_mask", "_labels", "_crop"):
        if stem.endswith(suffix):
            aliases.append(stem[: -len(suffix)])
    return aliases


def _build_mask_lookup(masks_dir: Path) -> dict[str, Path]:
    lookup: dict[str, Path] = {}
    for mask_path in sorted(masks_dir.glob("*.tif*")):
        for alias in _stem_aliases(mask_path.stem):
            lookup.setdefault(alias, mask_path)
    return lookup


def _find_crop(crops_dir: Path, sample: str) -> Path | None:
    for extension in (".tif", ".tiff", ".TIF", ".TIFF"):
        candidate = crops_dir / f"{sample}{extension}"
        if candidate.is_file():
            return candidate
    return None


def _find_mask(
    sample: str,
    lookup: dict[str, Path],
    masks_dir: Path,
) -> Path | None:
    for alias in _stem_aliases(sample):
        if alias in lookup:
            return lookup[alias]
    for suffix in ("", "_mask", "_labels"):
        for extension in (".tif", ".tiff", ".TIF", ".TIFF"):
            candidate = masks_dir / f"{sample}{suffix}{extension}"
            if candidate.is_file():
                return candidate
    return None


def _pair_crops_and_masks(
    crops_dir: Path,
    masks_dir: Path,
    *,
    existing: pd.DataFrame | None = None,
) -> list[tuple[Path, Path, str]]:
    if not crops_dir.is_dir():
        raise FileNotFoundError(f"Crops directory not found: {crops_dir}")
    if not masks_dir.is_dir():
        raise FileNotFoundError(f"Masks directory not found: {masks_dir}")

    lookup = _build_mask_lookup(masks_dir)
    crop_files = sorted(crops_dir.glob("*.tif*"))
    if not crop_files and existing is None:
        raise FileNotFoundError(f"No crop TIFFs found in {crops_dir}")

    samples: list[str]
    if existing is not None and "sample" in existing.columns:
        samples = existing["sample"].astype(str).tolist()
    else:
        samples = [crop_path.stem for crop_path in crop_files]

    pairs: list[tuple[Path, Path, str]] = []
    missing_crops: list[str] = []
    missing_masks: list[str] = []
    for sample in samples:
        crop_path = _find_crop(crops_dir, sample)
        if crop_path is None:
            missing_crops.append(sample)
            continue
        mask_path = _find_mask(sample, lookup, masks_dir)
        if mask_path is None:
            missing_masks.append(sample)
            continue
        pairs.append((crop_path, mask_path, sample))

    if not pairs:
        message = [
            f"No crop/mask pairs resolved under {crops_dir} and {masks_dir}.",
            f"Crop TIFF count: {len(crop_files)}",
            f"Mask TIFF count: {len(list(masks_dir.glob('*.tif*')))}",
        ]
        if missing_crops:
            message.append(
                "Missing crops (first 5): " + ", ".join(missing_crops[:5])
            )
        if missing_masks:
            message.append(
                "Missing masks (first 5): " + ", ".join(missing_masks[:5])
            )
        if crop_files and missing_masks:
            example_crop = crop_files[0].name
            example_masks = sorted(path.name for path in masks_dir.glob("*.tif*"))[:5]
            message.append(f"Example crop file: {example_crop}")
            message.append(
                "Example mask files: " + (", ".join(example_masks) or "(none)")
            )
        raise FileNotFoundError("\n".join(message))

    if missing_crops or missing_masks:
        print(
            f"Warning: skipped {len(missing_crops)} missing crops and "
            f"{len(missing_masks)} missing masks; continuing with {len(pairs)} pairs",
            file=sys.stderr,
        )

    return pairs


def measure_crop(crop_path: Path, mask_path: Path) -> dict[str, float]:
    from cellprofiler.modules.measureobjectintensity import MeasureObjectIntensity
    from cellprofiler.modules.measureobjectsizeshape import MeasureObjectSizeShape
    from cellprofiler_core.image import Image, ImageSetList
    from cellprofiler_core.measurement import Measurements
    from cellprofiler_core.object import ObjectSet, Objects
    from cellprofiler_core.pipeline import Pipeline
    from cellprofiler_core.preferences import set_headless
    from cellprofiler_core.workspace import Workspace

    image_hwc = _load_image(crop_path)
    labels = _labels_from_mask(np.asarray(tifffile.imread(mask_path)))
    if labels.shape != image_hwc.shape[:2]:
        raise ValueError(
            f"Crop/mask shape mismatch for {crop_path.name}: "
            f"{image_hwc.shape[:2]} vs {labels.shape}"
        )

    set_headless()
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
        if len(values) == 0:
            raise ValueError(f"No objects measured in {crop_path.name}")
        if len(values) != 1:
            raise ValueError(
                f"Expected one object in {crop_path.name}, got {len(values)}"
            )
        row[feature] = float(values[0])
    return row


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
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print pairing summary and exit without measuring",
    )
    args = parser.parse_args()

    analysis_root = args.analysis_root.expanduser()
    crops_dir = args.crops_dir or (analysis_root / "results" / "single_cell_crops")
    masks_dir = args.masks_dir or (analysis_root / "results" / "single_cell_masks")
    output_path = args.output or (analysis_root / "results" / "single_cell_full_features.csv")

    existing: pd.DataFrame | None = None
    existing_path = args.existing_features
    if existing_path is None:
        default_existing = analysis_root / "results" / "single_cell_features_and_clusters.csv"
        if default_existing.is_file():
            existing_path = default_existing
    if existing_path is not None:
        existing_path = existing_path.expanduser()
        if existing_path.is_file():
            existing = pd.read_csv(existing_path)
            print(f"Using existing feature table for pairing/metadata: {existing_path}")

    pairs = _pair_crops_and_masks(crops_dir, masks_dir, existing=existing)
    if args.limit is not None:
        pairs = pairs[: args.limit]

    print(f"Crops dir: {crops_dir}")
    print(f"Masks dir: {masks_dir}")
    print(f"Resolved {len(pairs)} crop/mask pairs")
    if pairs:
        crop_path, mask_path, sample = pairs[0]
        print(f"First pair: sample={sample}")
        print(f"  crop={crop_path.name}")
        print(f"  mask={mask_path.name}")

    if args.dry_run:
        return

    _require_cellprofiler()

    rows: list[dict[str, object]] = []
    feature_order: list[str] | None = None
    for index, (crop_path, mask_path, sample) in enumerate(pairs, start=1):
        try:
            measurements = measure_crop(crop_path, mask_path)
        except Exception as exc:
            raise RuntimeError(
                f"Failed measuring {crop_path.name} with mask {mask_path.name}: {exc}"
            ) from exc

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
                for column in (
                    "cluster",
                    "original_cell_id",
                    "compact_cell_id",
                    "crop_size_pixels",
                    "measurement_backend",
                ):
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
