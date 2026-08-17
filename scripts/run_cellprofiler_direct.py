"""Run pinned CellProfiler measurement modules on preloaded TIFF arrays."""

from __future__ import annotations

import argparse
import csv
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


def _scaled_plane(path: Path) -> np.ndarray:
    plane = np.asarray(tifffile.imread(path))
    if plane.ndim != 2:
        raise ValueError(f"Expected one grayscale plane at {path}, got {plane.shape}")
    if np.issubdtype(plane.dtype, np.integer):
        return plane.astype(np.float32) / float(np.iinfo(plane.dtype).max)
    result = plane.astype(np.float32)
    if np.nanmin(result) < 0 or np.nanmax(result) > 1:
        raise ValueError(f"Floating-point intensity image is outside [0, 1]: {path}")
    return result


def _measure(
    image_path: Path,
    labels_path: Path,
    image_name: str,
    object_name: str,
) -> dict[str, np.ndarray]:
    image_pixels = _scaled_plane(image_path)
    labels = np.asarray(tifffile.imread(labels_path))
    if labels.shape != image_pixels.shape:
        raise ValueError(
            f"Image/label shape mismatch: {image_pixels.shape} != {labels.shape}"
        )

    pipeline = Pipeline()
    image_set_list = ImageSetList()
    image_set = image_set_list.get_image_set(0)
    image_set.add(image_name, Image(image_pixels, dimensions=2))
    objects = Objects()
    objects.segmented = labels.astype(np.int32, copy=False)
    object_set = ObjectSet()
    object_set.add_objects(objects, object_name)
    measurements = Measurements()

    size_shape = MeasureObjectSizeShape()
    size_shape.set_module_num(1)
    size_shape.objects_list.value = object_name
    size_shape.calculate_zernikes.value = True
    size_shape.calculate_advanced.value = True
    pipeline.add_module(size_shape)

    intensity = MeasureObjectIntensity()
    intensity.set_module_num(2)
    intensity.images_list.value = image_name
    intensity.objects_list.value = object_name
    pipeline.add_module(intensity)

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

    return {
        feature: np.asarray(measurements.get_measurement(object_name, feature))
        for feature in measurements.get_feature_names(object_name)
    }


def run(
    load_data_path: Path,
    output_path: Path,
    image_name: str,
    object_name: str,
) -> None:
    set_headless()
    with load_data_path.open("r", encoding="utf-8-sig", newline="") as handle:
        inputs = list(csv.DictReader(handle))

    feature_order: list[str] | None = None
    output_rows: list[dict[str, object]] = []
    for image_number, input_row in enumerate(inputs, start=1):
        image_path = (
            Path(input_row[f"Image_PathName_{image_name}"])
            / input_row[f"Image_FileName_{image_name}"]
        )
        labels_path = (
            Path(input_row[f"Image_ObjectsPathName_{object_name}"])
            / input_row[f"Image_ObjectsFileName_{object_name}"]
        )
        measurements = _measure(image_path, labels_path, image_name, object_name)
        current_order = list(measurements)
        if feature_order is None:
            feature_order = current_order
        elif current_order != feature_order:
            raise ValueError("CellProfiler feature order changed between images")
        object_count = len(next(iter(measurements.values()), []))
        for object_index in range(object_count):
            row: dict[str, object] = {
                "ImageNumber": image_number,
                "ObjectNumber": object_index + 1,
            }
            row.update(
                {
                    feature: values[object_index]
                    for feature, values in measurements.items()
                }
            )
            output_rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["ImageNumber", "ObjectNumber", *(feature_order or [])]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--load-data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--image-name", required=True)
    parser.add_argument("--object-name", required=True)
    args = parser.parse_args()
    run(args.load_data, args.output, args.image_name, args.object_name)


if __name__ == "__main__":
    main()
