"""Measure full CellProfiler SizeShape+Intensity features on per-cell crop TIFFs."""

from __future__ import annotations

import argparse
import csv
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import tifffile
from cellprofiler.modules.measureobjectintensity import MeasureObjectIntensity
from cellprofiler.modules.measureobjectsizeshape import MeasureObjectSizeShape
from cellprofiler_core.image import Image, ImageSetList
from cellprofiler_core.measurement import Measurements
from cellprofiler_core.object import ObjectSet, Objects
from cellprofiler_core.pipeline import Pipeline
from cellprofiler_core.preferences import set_headless
from cellprofiler_core.workspace import Workspace
from scipy.ndimage import convolve

try:
    from skimage.morphology import skeletonize

    _HAS_SKIMAGE = True
except ImportError:
    _HAS_SKIMAGE = False

CELL_STEM = re.compile(r"^(?P<image_sample>.+)__cell_(?P<original_cell_id>\d+)$")
TIMEPOINT = re.compile(r"(?<![A-Za-z0-9])(ctrl|4h|24h|72h|7d)(?![A-Za-z0-9])", re.IGNORECASE)
CHECKPOINT_EVERY = 500
META_COLUMNS = ("cluster", "compact_cell_id")


def _channel_count(path: Path) -> int:
    with tifffile.TiffFile(path) as handle:
        array = handle.asarray(out="memmap")
    if array.ndim == 2:
        return 1
    if array.ndim == 3 and array.shape[0] <= 8 and array.shape[0] < array.shape[-1]:
        return int(array.shape[0])
    if array.ndim == 3:
        return int(array.shape[-1])
    raise ValueError(f"Unsupported crop shape {array.shape} at {path}")


def _max_channels(pairs: list[tuple[Path, Path, str, int]]) -> int:
    by_image: dict[str, Path] = {}
    for crop_path, _, image_sample, _ in pairs:
        by_image.setdefault(image_sample, crop_path)
    return max(_channel_count(path) for path in by_image.values())


def _load_image(path: Path, channel_count: int | None = None) -> np.ndarray:
    array = np.asarray(tifffile.imread(path))
    if array.ndim == 2:
        image = array[..., np.newaxis]
    elif array.ndim == 3 and array.shape[0] <= 8 and array.shape[0] < array.shape[-1]:
        image = np.moveaxis(array, 0, -1)
    elif array.ndim == 3:
        image = array
    else:
        raise ValueError(f"Unsupported crop shape {array.shape} at {path}")

    if channel_count is not None and image.shape[-1] < channel_count:
        pad_shape = (*image.shape[:2], channel_count - image.shape[-1])
        image = np.concatenate(
            [image, np.zeros(pad_shape, dtype=image.dtype)],
            axis=-1,
        )
    return image


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


def skeleton_features(cell_mask: np.ndarray) -> tuple[float, float]:
    if not _HAS_SKIMAGE:
        return float("nan"), float("nan")
    skeleton = skeletonize(cell_mask.astype(bool))
    neighbor_count = convolve(skeleton.astype(np.uint8), np.ones((3, 3), dtype=np.uint8))
    branch_points = skeleton & ((neighbor_count - skeleton.astype(bool)) >= 3)
    return float(skeleton.sum()), float(branch_points.sum())


def measure_crop(crop_path: Path, mask_path: Path, channel_count: int) -> dict[str, float]:
    image_hwc = _load_image(crop_path, channel_count=channel_count)
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
    intensity.images_list.value = ", ".join(channel_names)
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

    cell_mask = labels > 0
    skeleton_length, branch_points = skeleton_features(cell_mask)
    row["skeleton_length_pixels"] = skeleton_length
    row["nontrunk_branch_points"] = branch_points
    return row


def parse_crop_stem(stem: str) -> tuple[str, int]:
    match = CELL_STEM.match(stem)
    if match is None:
        raise ValueError(f"Unexpected crop filename stem: {stem}")
    return match.group("image_sample"), int(match.group("original_cell_id"))


def _pair_crops_and_masks(
    crops_dir: Path, masks_dir: Path
) -> list[tuple[Path, Path, str, int]]:
    pairs: list[tuple[Path, Path, str, int]] = []
    for crop_path in sorted(crops_dir.glob("*.tif*")):
        image_sample, original_cell_id = parse_crop_stem(crop_path.stem)
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
        pairs.append((crop_path, mask_path, image_sample, original_cell_id))
    if not pairs:
        raise FileNotFoundError(f"No crop TIFFs found in {crops_dir}")
    return pairs


def _parse_time(image_sample: str) -> str | None:
    match = TIMEPOINT.search(image_sample)
    return match.group(1).lower() if match else None


def _read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return [], []
        return list(reader.fieldnames), list(reader)


def _metadata_lookup(existing_path: Path | None) -> dict[tuple[str, int], dict[str, str]]:
    if existing_path is None or not existing_path.is_file():
        return {}
    _, rows = _read_csv_rows(existing_path)
    lookup: dict[tuple[str, int], dict[str, str]] = {}
    for row in rows:
        if "sample" not in row or "original_cell_id" not in row:
            continue
        key = (row["sample"], int(row["original_cell_id"]))
        lookup[key] = row
    return lookup


def _done_crop_paths(output_path: Path) -> set[str]:
    if not output_path.is_file():
        return set()
    _, rows = _read_csv_rows(output_path)
    return {row["crop_path"] for row in rows if row.get("crop_path")}


def _feature_columns(rows: list[dict[str, object]]) -> list[str]:
    area_shape: list[str] = []
    intensity: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key in seen:
                continue
            if key.startswith("AreaShape_"):
                area_shape.append(key)
                seen.add(key)
            elif key.startswith("Intensity_"):
                intensity.append(key)
                seen.add(key)
    return sorted(area_shape) + sorted(intensity)


def _row_fieldnames(rows: list[dict[str, object]]) -> list[str]:
    leading = [
        "sample",
        "original_cell_id",
        "compact_cell_id",
        "Time",
        "cluster",
        "measurement_backend",
        "crop_path",
        "mask_path",
    ]
    features = _feature_columns(rows)
    trailing = ["skeleton_length_pixels", "nontrunk_branch_points"]
    fieldnames = [column for column in leading if any(column in row for row in rows)]
    for column in features + trailing:
        if column not in fieldnames:
            fieldnames.append(column)
    return fieldnames


def _write_table(rows: list[dict[str, object]], output_path: Path) -> None:
    if not rows:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = _row_fieldnames(rows)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _measure_worker(
    task: tuple[Path, Path, str, int, int],
) -> dict[str, object]:
    crop_path, mask_path, image_sample, original_cell_id, channel_count = task
    measurements = measure_crop(crop_path, mask_path, channel_count)
    return {
        "sample": image_sample,
        "original_cell_id": original_cell_id,
        "crop_path": str(crop_path),
        "mask_path": str(mask_path),
        "Time": _parse_time(image_sample) or "unknown",
        "measurement_backend": "cellprofiler",
        **measurements,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis"),
    )
    parser.add_argument("--crops-dir", type=Path, default=None)
    parser.add_argument("--masks-dir", type=Path, default=None)
    parser.add_argument("--existing-features", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip cells already present in the output CSV (by crop_path).",
    )
    args = parser.parse_args()

    analysis_root = args.analysis_root.expanduser()
    crops_dir = args.crops_dir or (analysis_root / "results" / "single_cell_crops")
    masks_dir = args.masks_dir or (analysis_root / "results" / "single_cell_masks")
    output_path = args.output or (analysis_root / "results" / "single_cell_full_features.csv")

    set_headless()
    pairs = _pair_crops_and_masks(crops_dir, masks_dir)
    channel_count = _max_channels(pairs)
    metadata = _metadata_lookup(
        args.existing_features.expanduser() if args.existing_features else None
    )

    if args.skip_existing:
        done_paths = _done_crop_paths(output_path)
        pairs = [pair for pair in pairs if str(pair[0]) not in done_paths]

    if args.limit is not None:
        pairs = pairs[: args.limit]

    resumed_rows: list[dict[str, object]] = []
    if args.skip_existing and output_path.is_file():
        _, resumed_rows_raw = _read_csv_rows(output_path)
        resumed_rows = [dict(row) for row in resumed_rows_raw]

    rows: list[dict[str, object]] = list(resumed_rows)
    total = len(pairs)
    workers = max(1, args.workers)

    if total == 0:
        print(f"No new cells to measure; output already has {len(rows)} rows at {output_path}")
        return

    tasks = [
        (crop_path, mask_path, image_sample, original_cell_id, channel_count)
        for crop_path, mask_path, image_sample, original_cell_id in pairs
    ]
    print(f"Measuring {total} cells with {workers} worker(s) across {channel_count} channel slot(s)")
    if not _HAS_SKIMAGE:
        print("WARNING: scikit-image unavailable; skeleton columns will be NaN")

    with Pool(processes=workers) as pool:
        for index, row in enumerate(pool.imap_unordered(_measure_worker, tasks, chunksize=8), start=1):
            meta = metadata.get((str(row["sample"]), int(row["original_cell_id"])))
            if meta is not None:
                for column in META_COLUMNS:
                    if meta.get(column):
                        row[column] = meta[column]

            rows.append(row)
            if index % 250 == 0 or index == total:
                print(f"Measured {index}/{total} new cells")
            if index % CHECKPOINT_EVERY == 0 or index == total:
                _write_table(rows, output_path)

    _write_table(rows, output_path)
    feature_count = len(_feature_columns(rows)) + 2
    print(
        f"Wrote {len(rows)} rows x {feature_count} features "
        f"(CellProfiler + skeleton) to {output_path}"
    )


if __name__ == "__main__":
    main()
