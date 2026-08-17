"""Exact CellProfiler-backed measurement execution and output assembly."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import tifffile

from .images import load_hwc
from .inventory import ImageRecord
from .schema import MetricSpec, SchemaError, inspect_cppipe, read_csv_header


class CellProfilerError(RuntimeError):
    """Raised when exact CellProfiler measurement cannot be guaranteed."""


def cellprofiler_version(executable: str) -> str:
    """Return the installed CellProfiler version string."""
    resolved = shutil.which(executable)
    if resolved is None:
        raise CellProfilerError(f"CellProfiler executable was not found: {executable}")
    completed = subprocess.run(
        [resolved, "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise CellProfilerError(
            f"Could not query CellProfiler version: {completed.stderr.strip()}"
        )
    match = re.search(r"(\d+\.\d+(?:\.\d+)*)", completed.stdout + completed.stderr)
    if not match:
        raise CellProfilerError(
            f"Could not parse CellProfiler version from: {completed.stdout.strip()}"
        )
    return match.group(1)


def prepare_load_data(
    records: list[ImageRecord],
    mask_paths: dict[str, Path],
    destination: Path,
    image_name: str,
    object_name: str,
    image_channel_index: int | None = None,
) -> pd.DataFrame:
    """Write a deterministic LoadData CSV that injects Cellpose labels as objects."""
    rows: list[dict[str, Any]] = []
    input_dir = destination.parent / "inputs"
    for image_number, record in enumerate(records, start=1):
        image_path = Path(record.path).resolve()
        if image_channel_index is not None:
            if image_channel_index < 0 or image_channel_index >= record.channel_count:
                raise CellProfilerError(
                    f"Invalid CellProfiler channel index {image_channel_index} "
                    f"for {record.image_id}"
                )
            plane = load_hwc(record)[..., image_channel_index]
            input_dir.mkdir(parents=True, exist_ok=True)
            image_path = (input_dir / f"{record.image_id}_{image_name}.tif").resolve()
            tifffile.imwrite(
                image_path,
                plane,
                metadata={"axes": "YX"},
                photometric="minisblack",
            )
        mask_path = mask_paths[record.image_id].resolve()
        rows.append(
            {
                f"Image_FileName_{image_name}": image_path.name,
                f"Image_PathName_{image_name}": str(image_path.parent),
                f"Image_ObjectsFileName_{object_name}": mask_path.name,
                f"Image_ObjectsPathName_{object_name}": str(mask_path.parent),
                "Metadata_ImageID": record.image_id,
                "Metadata_SourceImageNumber": image_number,
            }
        )
    table = pd.DataFrame.from_records(rows)
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, index=False)
    return table


def _required_modules(metrics: Iterable[MetricSpec]) -> set[str]:
    return {
        metric.module
        for metric in metrics
        if metric.role == "feature"
        and metric.module
        not in {"Unmapped", "IdentifyObjects/MeasureObjectIntensity", "IdentifyObjects"}
    }


def validate_measurement_pipeline(
    pipeline_path: Path,
    metrics: list[MetricSpec],
) -> None:
    """Ensure the selected pipeline can load external objects and measure the schema."""
    if not pipeline_path.is_file():
        raise CellProfilerError(f"Measurement pipeline does not exist: {pipeline_path}")
    spec = inspect_cppipe(pipeline_path)
    if "LoadData" not in spec.modules:
        raise CellProfilerError(
            "Exact Cellpose-object injection requires a measurement pipeline with LoadData. "
            "Export a copy of the original pipeline that loads the supplied image and "
            "Image_ObjectsFileName object columns while preserving all measurement settings."
        )
    missing = sorted(_required_modules(metrics) - set(spec.modules))
    if missing:
        raise CellProfilerError(
            f"Measurement pipeline is missing required modules: {', '.join(missing)}"
        )
    if not any(module.startswith("ExportTo") for module in spec.modules):
        raise CellProfilerError("Measurement pipeline must contain an ExportTo* module")


def _find_matching_export(output_dir: Path, requested_columns: list[str]) -> Path:
    candidates: list[tuple[int, Path]] = []
    feature_columns = [
        column
        for column in requested_columns
        if column not in {"ImageNumber", "ObjectNumber"} and not column.startswith("Metadata_")
    ]
    for path in output_dir.glob("*.csv"):
        try:
            header = set(read_csv_header(path))
        except (SchemaError, UnicodeError):
            continue
        score = sum(column in header for column in feature_columns)
        candidates.append((score, path))
    if not candidates:
        raise CellProfilerError(f"CellProfiler produced no readable CSV files in {output_dir}")
    score, path = max(candidates, key=lambda item: (item[0], item[1].name))
    if score != len(feature_columns):
        missing = sorted(set(feature_columns) - set(read_csv_header(path)))
        raise CellProfilerError(
            f"No CellProfiler export matched the requested schema. Best file {path.name} "
            f"is missing: {missing}"
        )
    return path


def _metadata_value(record: ImageRecord, column: str) -> Any:
    base = re.sub(r"\.\d+$", "", column)
    filename = Path(record.path).name
    if base == "Metadata_Time":
        return filename.split("_", maxsplit=1)[0]
    if base == "Metadata_PhotoID":
        match = re.search(r"_(\d+)\.(?:ome\.)?tiff?$", filename, flags=re.IGNORECASE)
        return match.group(1) if match else pd.NA
    if base == "Metadata_FileLocation":
        return str(Path(record.path).parent)
    if base == "Metadata_SizeY":
        return record.height
    if base == "Metadata_SizeX":
        return record.width
    if base == "Metadata_SizeC":
        return record.channel_count
    return pd.NA


def _add_requested_metadata(
    table: pd.DataFrame,
    requested_columns: list[str],
    records: list[ImageRecord],
) -> pd.DataFrame:
    if "ImageNumber" not in table:
        raise CellProfilerError("CellProfiler export is missing ImageNumber")
    record_lookup = {
        image_number: record for image_number, record in enumerate(records, start=1)
    }
    output = table.copy()
    for column in requested_columns:
        if column in output.columns:
            continue
        if not column.startswith("Metadata_"):
            continue
        output[column] = [
            _metadata_value(record_lookup.get(int(image_number)), column)
            if int(image_number) in record_lookup
            else pd.NA
            for image_number in output["ImageNumber"]
        ]
    return output


def _cellprofiler_package_version(python_executable: str) -> str:
    completed = subprocess.run(
        [
            python_executable,
            "-c",
            "import importlib.metadata; print(importlib.metadata.version('cellprofiler'))",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise CellProfilerError(
            "Could not query CellProfiler package version: "
            + completed.stderr.strip()
        )
    return completed.stdout.strip()


def _run_cellprofiler_direct(
    config: dict[str, Any],
    records: list[ImageRecord],
    mask_paths: dict[str, Path],
    metrics: list[MetricSpec],
    output_root: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cp_config = config["cellprofiler"]
    pipeline_path = Path(config["paths"]["cellprofiler_measurement_pipeline"]).expanduser()
    validate_measurement_pipeline(pipeline_path, metrics)
    python_executable = str(cp_config["python_executable"])
    direct_script = Path(cp_config["direct_script"]).expanduser()
    if not Path(python_executable).is_file():
        raise CellProfilerError(
            f"CellProfiler Python executable does not exist: {python_executable}"
        )
    if not direct_script.is_file():
        raise CellProfilerError(f"Direct measurement script does not exist: {direct_script}")
    installed_version = _cellprofiler_package_version(python_executable)
    expected_version = cp_config.get("expected_version")
    if installed_version != str(expected_version):
        raise CellProfilerError(
            f"CellProfiler version mismatch: expected {expected_version}, "
            f"found {installed_version}"
        )

    cp_root = output_root / "cellprofiler"
    cp_output = cp_root / "output"
    cp_output.mkdir(parents=True, exist_ok=True)
    load_data_path = cp_root / "load_data.csv"
    prepare_load_data(
        records,
        mask_paths,
        load_data_path,
        image_name=str(cp_config["image_name"]),
        object_name=str(cp_config["object_name"]),
        image_channel_index=int(cp_config["image_channel_index"]),
    )
    export_path = cp_output / "AstroResultsCellpose_masks.csv"
    command = [
        python_executable,
        str(direct_script),
        "--load-data",
        str(load_data_path),
        "--output",
        str(export_path),
        "--image-name",
        str(cp_config["image_name"]),
        "--object-name",
        str(cp_config["object_name"]),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    log_path = cp_root / "cellprofiler.log"
    log_path.write_text(
        f"COMMAND: {command!r}\n\nSTDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise CellProfilerError(
            f"Direct CellProfiler measurement failed with exit code "
            f"{completed.returncode}; see {log_path}"
        )
    table = pd.read_csv(export_path)
    requested_columns = [metric.column for metric in metrics]
    required_features = [metric.column for metric in metrics if metric.role == "feature"]
    missing_features = [
        column for column in required_features if column not in table.columns
    ]
    if missing_features:
        raise CellProfilerError(
            f"CellProfiler direct output is missing requested features: {missing_features}"
        )
    table = _add_requested_metadata(table, requested_columns, records)
    provenance = {
        "backend": "direct",
        "cellprofiler_version": installed_version,
        "pipeline_path": str(pipeline_path),
        "pipeline": inspect_cppipe(pipeline_path).__dict__,
        "load_data_path": str(load_data_path),
        "export_path": str(export_path),
        "command": command,
    }
    (cp_root / "run_manifest.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True, default=list),
        encoding="utf-8",
    )
    return table[requested_columns].copy(), provenance


def run_cellprofiler(
    config: dict[str, Any],
    records: list[ImageRecord],
    mask_paths: dict[str, Path],
    metrics: list[MetricSpec],
    output_root: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run version-pinned CellProfiler headlessly on Cellpose label objects."""
    cp_config = config["cellprofiler"]
    if cp_config.get("backend", "cli") == "direct":
        return _run_cellprofiler_direct(
            config, records, mask_paths, metrics, output_root
        )
    pipeline_value = config["paths"].get("cellprofiler_measurement_pipeline")
    if not pipeline_value:
        raise CellProfilerError(
            "Set paths.cellprofiler_measurement_pipeline after inventory identifies "
            "the version-matched pipeline that loads Cellpose labels."
        )
    pipeline_path = Path(pipeline_value).expanduser()
    validate_measurement_pipeline(pipeline_path, metrics)
    executable = str(cp_config["executable"])
    installed_version = cellprofiler_version(executable)
    expected_version = cp_config.get("expected_version")
    if not expected_version and not cp_config.get("allow_unpinned_version", False):
        raise CellProfilerError(
            "cellprofiler.expected_version must be set for exact measurements"
        )
    if expected_version and installed_version != str(expected_version):
        raise CellProfilerError(
            f"CellProfiler version mismatch: expected {expected_version}, found {installed_version}"
        )
    cp_root = output_root / "cellprofiler"
    cp_output = cp_root / "output"
    cp_output.mkdir(parents=True, exist_ok=True)
    load_data_path = cp_root / "load_data.csv"
    prepare_load_data(
        records,
        mask_paths,
        load_data_path,
        image_name=str(cp_config["image_name"]),
        object_name=str(cp_config["object_name"]),
        image_channel_index=int(cp_config["image_channel_index"]),
    )
    resolved = shutil.which(executable)
    assert resolved is not None  # Checked by cellprofiler_version.
    command = [
        resolved,
        "--run-headless",
        "--run",
        "--pipeline",
        str(pipeline_path),
        "--data-file",
        str(load_data_path),
        "--output-directory",
        str(cp_output),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    log_path = cp_root / "cellprofiler.log"
    log_path.write_text(
        f"COMMAND: {command!r}\n\nSTDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise CellProfilerError(
            f"CellProfiler failed with exit code {completed.returncode}; see {log_path}"
        )
    requested_columns = [metric.column for metric in metrics]
    export_path = _find_matching_export(cp_output, requested_columns)
    table = pd.read_csv(export_path)
    required_features = [metric.column for metric in metrics if metric.role == "feature"]
    missing_features = [
        column for column in required_features if column not in table.columns
    ]
    if missing_features:
        raise CellProfilerError(
            f"CellProfiler export is missing requested features: {missing_features}"
        )
    table = _add_requested_metadata(table, requested_columns, records)
    missing = [column for column in requested_columns if column not in table.columns]
    if missing:
        raise CellProfilerError(f"CellProfiler export is missing requested columns: {missing}")
    provenance = {
        "cellprofiler_version": installed_version,
        "pipeline_path": str(pipeline_path),
        "pipeline": inspect_cppipe(pipeline_path).__dict__,
        "load_data_path": str(load_data_path),
        "export_path": str(export_path),
        "command": command,
    }
    (cp_root / "run_manifest.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True, default=list),
        encoding="utf-8",
    )
    return table[requested_columns].copy(), provenance


def measurement_family(table: pd.DataFrame, category: str) -> pd.DataFrame:
    """Select one CellProfiler measurement family without recomputing values."""
    prefix = f"{category}_"
    columns = [column for column in table.columns if column.startswith(prefix)]
    return table[columns].copy()


def size_shape_features(table: pd.DataFrame) -> pd.DataFrame:
    return measurement_family(table, "AreaShape")


def intensity_features(table: pd.DataFrame) -> pd.DataFrame:
    return measurement_family(table, "Intensity")


def texture_features(table: pd.DataFrame) -> pd.DataFrame:
    return measurement_family(table, "Texture")


def neighbor_features(table: pd.DataFrame) -> pd.DataFrame:
    return measurement_family(table, "Neighbors")


def assemble_cell_table(
    measurements: pd.DataFrame,
    nodes: pd.DataFrame,
    records: list[ImageRecord],
    requested_columns: list[str],
    output_root: Path,
) -> pd.DataFrame:
    """Join exact measurements to global coordinates without losing cells."""
    if "ImageNumber" not in measurements or "ObjectNumber" not in measurements:
        raise CellProfilerError("CellProfiler object export requires ImageNumber and ObjectNumber")
    image_lookup = {
        image_number: record.image_id for image_number, record in enumerate(records, start=1)
    }
    measured = measurements.copy()
    measured["image_id"] = measured["ImageNumber"].map(image_lookup)
    if measured["image_id"].isna().any():
        raise CellProfilerError("CellProfiler returned an unknown ImageNumber")
    measured["cell_id"] = measured["ObjectNumber"].astype(int)
    merged = nodes.merge(
        measured,
        on=["image_id", "cell_id"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if not (merged["_merge"] == "both").all():
        counts = merged["_merge"].value_counts().to_dict()
        raise CellProfilerError(
            "Every Cellpose object must have exactly one CellProfiler vector; "
            f"join result={counts}"
        )
    merged = merged.drop(columns="_merge")
    identity = ["image_id", "cell_id", "global_cell_id"]
    metric_order = [column for column in requested_columns if column not in identity]
    spatial_order = [
        column for column in nodes.columns if column not in {"image_id", "cell_id", "global_cell_id"}
    ]
    ordered = merged[identity + metric_order + spatial_order]
    table_dir = output_root / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    ordered.to_csv(table_dir / "cells.csv", index=False)
    ordered.to_parquet(table_dir / "cells.parquet", index=False)
    return ordered


def assert_numeric_parity(
    reference: pd.DataFrame,
    candidate: pd.DataFrame,
    key_columns: list[str],
    rtol: float = 1.0e-7,
    atol: float = 1.0e-9,
) -> None:
    """Assert keyed schema and numerical parity against a golden CP export."""
    if list(reference.columns) != list(candidate.columns):
        raise AssertionError("CellProfiler parity schema/order mismatch")
    left = reference.sort_values(key_columns).reset_index(drop=True)
    right = candidate.sort_values(key_columns).reset_index(drop=True)
    if len(left) != len(right):
        raise AssertionError(f"CellProfiler parity row mismatch: {len(left)} != {len(right)}")
    for column in left.columns:
        if pd.api.types.is_numeric_dtype(left[column]):
            np.testing.assert_allclose(
                left[column].to_numpy(),
                right[column].to_numpy(),
                rtol=rtol,
                atol=atol,
                equal_nan=True,
                err_msg=f"CellProfiler parity failed for {column}",
            )
        elif not left[column].equals(right[column]):
            raise AssertionError(f"CellProfiler parity failed for {column}")
