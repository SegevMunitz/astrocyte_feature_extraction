from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from astrocyte_feature_extraction.cellprofiler_features import (
    assemble_cell_table,
    assert_numeric_parity,
    prepare_load_data,
)

from conftest import write_test_image


def test_load_data_injects_integer_objects(tmp_path: Path) -> None:
    record = write_test_image(tmp_path / "image.tif")
    mask = tmp_path / "labels.tif"
    mask.write_bytes(b"placeholder")
    destination = tmp_path / "load_data.csv"
    table = prepare_load_data(
        [record],
        {record.image_id: mask},
        destination,
        image_name="Input",
        object_name="Astrocytes",
        image_channel_index=0,
    )
    assert "Image_ObjectsFileName_Astrocytes" in table
    assert table.iloc[0]["Metadata_ImageID"] == "image_a"
    assert destination.is_file()
    plane_path = (
        Path(table.iloc[0]["Image_PathName_Input"])
        / table.iloc[0]["Image_FileName_Input"]
    )
    assert tifffile.imread(plane_path).shape == (record.height, record.width)


def test_assemble_preserves_schema_and_spatial_columns(tmp_path: Path) -> None:
    record = write_test_image(tmp_path / "image.tif")
    measurements = pd.DataFrame(
        {
            "ImageNumber": [1, 1],
            "ObjectNumber": [1, 2],
            "AreaShape_Area": [10.0, 20.0],
        }
    )
    nodes = pd.DataFrame(
        {
            "image_id": ["image_a", "image_a"],
            "cell_id": [1, 2],
            "global_cell_id": ["image_a:1", "image_a:2"],
            "centroid_x_px": [3.0, 8.0],
        }
    )
    table = assemble_cell_table(
        measurements,
        nodes,
        [record],
        ["ImageNumber", "ObjectNumber", "AreaShape_Area"],
        tmp_path,
    )
    assert list(table.columns) == [
        "image_id",
        "cell_id",
        "global_cell_id",
        "ImageNumber",
        "ObjectNumber",
        "AreaShape_Area",
        "centroid_x_px",
    ]
    assert (tmp_path / "tables" / "cells.parquet").is_file()


def test_numeric_parity_is_keyed() -> None:
    reference = pd.DataFrame(
        {"ImageNumber": [1, 1], "ObjectNumber": [1, 2], "AreaShape_Area": [4.0, 8.0]}
    )
    candidate = reference.iloc[::-1].reset_index(drop=True)
    assert_numeric_parity(
        reference, candidate, key_columns=["ImageNumber", "ObjectNumber"]
    )
    changed = candidate.copy()
    changed.loc[0, "AreaShape_Area"] += 1
    try:
        assert_numeric_parity(
            reference, changed, key_columns=["ImageNumber", "ObjectNumber"]
        )
    except AssertionError:
        pass
    else:
        raise AssertionError("Expected golden parity mismatch")
