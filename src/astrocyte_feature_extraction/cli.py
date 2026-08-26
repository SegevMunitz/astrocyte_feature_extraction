"""Command-line orchestration for the complete astrocyte pipeline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .cellprofiler_features import (
    CellProfilerError,
    assemble_cell_table,
    run_cellprofiler,
)
from .config import ConfigError, load_config, output_path
from .crop_measure import CropMeasureError, run_crop_measure, run_predict_rna
from .crops import export_all_crops
from .inventory import ImageRecord, discover_images, write_inventory
from .provenance import write_run_manifest
from .schema import (
    MetricSpec,
    PipelineSpec,
    build_metric_manifest,
    discover_pipelines,
    read_csv_header,
    write_schema_manifest,
)
from .segmentation import segment_all
from .spatial import build_spatial_tables

def _context(
    config: dict[str, Any],
) -> tuple[list[ImageRecord], list[MetricSpec], list[PipelineSpec], Path]:
    output_root = output_path(config)
    records = discover_images(config)
    metrics = build_metric_manifest(
        read_csv_header(config["paths"]["cellprofiler_results_csv"])
    )
    pipelines = discover_pipelines(config["paths"]["cellprofiler_pipeline_dir"])
    return records, metrics, pipelines, output_root


def _write_inventory_outputs(
    records: list[ImageRecord],
    metrics: list[MetricSpec],
    pipelines: list[PipelineSpec],
    output_root: Path,
) -> None:
    write_inventory(records, output_root / "manifests" / "images.json")
    write_schema_manifest(
        metrics, pipelines, output_root / "manifests" / "cellprofiler_schema.json"
    )


def _mask_paths(records: list[ImageRecord], output_root: Path) -> dict[str, Path]:
    paths = {
        record.image_id: output_root / "masks" / f"{record.image_id}_labels.tif"
        for record in records
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing segmentation masks; run the segment stage first: " + ", ".join(missing)
        )
    return paths


def _select_measurement_pipeline(
    config: dict[str, Any],
    metrics: list[MetricSpec],
    pipelines: list[PipelineSpec],
) -> None:
    if config["paths"].get("cellprofiler_measurement_pipeline"):
        return
    required = {
        metric.module
        for metric in metrics
        if metric.role == "feature"
        and metric.module
        not in {"Unmapped", "IdentifyObjects", "IdentifyObjects/MeasureObjectIntensity"}
    }
    candidates = [
        pipeline
        for pipeline in pipelines
        if "LoadData" in pipeline.modules
        and any(module.startswith("ExportTo") for module in pipeline.modules)
        and required.issubset(set(pipeline.modules))
    ]
    if len(candidates) != 1:
        paths = [pipeline.path for pipeline in candidates]
        raise CellProfilerError(
            "Could not uniquely select a Cellpose-object measurement pipeline. "
            f"Set paths.cellprofiler_measurement_pipeline explicitly; candidates={paths}"
        )
    config["paths"]["cellprofiler_measurement_pipeline"] = candidates[0].path


def run_inventory(config_path: Path) -> None:
    config = load_config(config_path)
    records, metrics, pipelines, output_root = _context(config)
    _write_inventory_outputs(records, metrics, pipelines, output_root)
    write_run_manifest(config, records, output_root, config_path)
    print(
        json.dumps(
            {
                "images": len(records),
                "requested_columns": len(metrics),
                "pipelines": len(pipelines),
                "output": str(output_root),
            },
            indent=2,
        )
    )


def run_segment(config_path: Path) -> None:
    config = load_config(config_path)
    records, metrics, pipelines, output_root = _context(config)
    _write_inventory_outputs(records, metrics, pipelines, output_root)
    masks = segment_all(records, config, output_root)
    print(f"Segmented {len(masks)} images into {output_root / 'masks'}")


def run_crops(config_path: Path) -> None:
    config = load_config(config_path)
    records, _, _, output_root = _context(config)
    table = export_all_crops(
        records, _mask_paths(records, output_root), output_root, config["crops"]
    )
    print(f"Exported {len(table)} isolated cells")


def run_spatial(config_path: Path) -> None:
    config = load_config(config_path)
    records, _, _, output_root = _context(config)
    nodes, edges = build_spatial_tables(
        records,
        _mask_paths(records, output_root),
        output_root,
        config["spatial"],
        config["qc"],
    )
    print(f"Exported {len(nodes)} nodes and {len(edges)} directed edges")


def run_measure(config_path: Path) -> None:
    config = load_config(config_path)
    records, metrics, pipelines, output_root = _context(config)
    masks = _mask_paths(records, output_root)
    _select_measurement_pipeline(config, metrics, pipelines)
    nodes, _ = build_spatial_tables(
        records, masks, output_root, config["spatial"], config["qc"]
    )
    measurements, _ = run_cellprofiler(
        config, records, masks, metrics, output_root
    )
    table = assemble_cell_table(
        measurements,
        nodes,
        records,
        [metric.column for metric in metrics],
        output_root,
    )
    print(f"Exported {len(table)} exact per-cell feature vectors")


def run_all(config_path: Path) -> None:
    config = load_config(config_path)
    records, metrics, pipelines, output_root = _context(config)
    _write_inventory_outputs(records, metrics, pipelines, output_root)
    write_run_manifest(config, records, output_root, config_path)
    masks = segment_all(records, config, output_root)
    export_all_crops(records, masks, output_root, config["crops"])
    morph_csv = run_crop_measure(config)
    rna_dir = run_predict_rna(config)
    print(
        f"Complete: {len(records)} images → morph {morph_csv.name} → RNA {rna_dir}"
    )


def run_crop_measure_stage(config_path: Path) -> None:
    config = load_config(config_path)
    morph_csv = run_crop_measure(config)
    print(f"Wrote crop-level morph features: {morph_csv}")


def run_predict_rna_stage(config_path: Path) -> None:
    config = load_config(config_path)
    out_dir = run_predict_rna(config)
    print(f"Wrote RNA predictions: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="astrocyte-pipeline",
        description=(
            "Images → Cellpose → crops → crop CellProfiler/skeleton → RNA predict."
        ),
    )
    parser.add_argument(
        "stage",
        choices=(
            "inventory",
            "segment",
            "crops",
            "spatial",
            "measure",
            "crop_measure",
            "predict_rna",
            "run",
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/elsc_unified.yaml"),
        help="YAML pipeline configuration",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    functions = {
        "inventory": run_inventory,
        "segment": run_segment,
        "crops": run_crops,
        "spatial": run_spatial,
        "measure": run_measure,
        "crop_measure": run_crop_measure_stage,
        "predict_rna": run_predict_rna_stage,
        "run": run_all,
    }
    try:
        functions[args.stage](args.config)
    except (ConfigError, CropMeasureError, RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
