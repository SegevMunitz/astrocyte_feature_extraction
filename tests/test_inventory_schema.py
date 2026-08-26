from __future__ import annotations

from pathlib import Path

from astrocyte_feature_extraction.inventory import inspect_image
from astrocyte_feature_extraction.schema import (
    build_metric_manifest,
    discover_pipelines,
    read_csv_header,
    validate_pipeline_modules,
)

from conftest import write_test_image


def test_inventory_reads_channel_axis(tmp_path: Path) -> None:
    image_path = tmp_path / "field.ome.tif"
    write_test_image(image_path)
    record = inspect_image(image_path, tmp_path, "auto", expected_channels=3)
    assert record.channel_axis == 0
    assert record.channel_count == 3
    assert (record.height, record.width) == (12, 16)


def test_discover_skips_incompatible_channels(tmp_path: Path) -> None:
    import numpy as np
    import tifffile
    from astrocyte_feature_extraction.inventory import discover_images

    good = tmp_path / "good.tif"
    write_test_image(good)
    bad = tmp_path / "bad.tif"
    tifffile.imwrite(
        bad,
        np.zeros((2, 8, 8), dtype=np.uint16),
        metadata={"axes": "CYX"},
        photometric="minisblack",
    )
    out = tmp_path / "out"
    out.mkdir()
    config = {
        "paths": {"images_dir": str(tmp_path), "output_dir": str(out)},
        "images": {
            "extensions": [".tif"],
            "recursive": False,
            "channel_axis": "auto",
            "expected_channels": 3,
            "skip_incompatible": True,
            "include_regex": None,
            "max_images": None,
        },
    }
    records = discover_images(config)
    assert len(records) == 1
    skipped = out / "manifests" / "skipped_images.json"
    assert skipped.is_file()


def test_schema_maps_columns_and_pipeline_modules(tmp_path: Path) -> None:
    csv_path = tmp_path / "objects.csv"
    csv_path.write_text(
        "ImageNumber,ObjectNumber,Metadata_Time,Metadata_Time,AreaShape_Area,"
        "Intensity_MeanIntensity_GFAP,Neighbors_NumberOfNeighbors_Adjacent,\n"
        "1,1,ctrl,ctrl,12,0.5,2,\n",
        encoding="utf-8",
    )
    columns = read_csv_header(csv_path)
    metrics = build_metric_manifest(columns)
    assert columns[2:4] == ["Metadata_Time", "Metadata_Time.1"]
    assert [metric.module for metric in metrics[-3:]] == [
        "MeasureObjectSizeShape",
        "MeasureObjectIntensity",
        "MeasureObjectNeighbors",
    ]

    pipeline_dir = tmp_path / "pipelines"
    pipeline_dir.mkdir()
    (pipeline_dir / "measure.cppipe").write_text(
        "CellProfiler Pipeline: http://www.cellprofiler.org\n"
        "Version:5\n"
        "DateRevision:424\n"
        "GitHash:abc\n"
        "ModuleCount:5\n"
        "LoadData:[module_num:1|variable_revision_number:6|]\n"
        "MeasureObjectSizeShape:[module_num:2|variable_revision_number:3|]\n"
        "MeasureObjectIntensity:[module_num:3|variable_revision_number:4|]\n"
        "MeasureObjectNeighbors:[module_num:4|variable_revision_number:3|]\n"
        "ExportToSpreadsheet:[module_num:5|variable_revision_number:13|]\n",
        encoding="utf-8",
    )
    pipelines = discover_pipelines(pipeline_dir)
    candidates = validate_pipeline_modules(metrics, pipelines)
    assert candidates["MeasureObjectSizeShape"] == [str(pipeline_dir / "measure.cppipe")]
