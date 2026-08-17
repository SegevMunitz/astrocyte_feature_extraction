from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from astrocyte_feature_extraction.crops import article_preprocess, export_cell_crops
from astrocyte_feature_extraction.segmentation import segment_record

from conftest import write_test_image


class FakeModel:
    input_image: np.ndarray | None = None

    def eval(
        self,
        image: np.ndarray,
        channels: Any = None,
        channel_axis: int = -1,
        diameter: float | None = None,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0.0,
        normalize: bool = True,
        min_size: int = 15,
    ) -> tuple[np.ndarray, None, None, None]:
        self.input_image = image
        labels = np.zeros(image.shape[:2], dtype=np.int32)
        labels[2:5, 3:7] = 4
        labels[7:10, 10:14] = 9
        return labels, None, None, None


def test_segment_and_export_cell_only_crops(tmp_path: Path, monkeypatch: Any) -> None:
    image_path = tmp_path / "image.tif"
    record = write_test_image(image_path)
    checkpoint = tmp_path / "weights"
    checkpoint.write_bytes(b"test checkpoint")
    config = {
        "paths": {"checkpoint": str(checkpoint)},
        "cellpose": {
            "expected_version": None,
            "gpu": False,
            "input_channel_indices": [1, 2, 0],
            "normalize": True,
            "diameter": None,
            "flow_threshold": 0.4,
            "cellprob_threshold": 0.0,
            "min_size": 1,
        },
    }
    monkeypatch.setattr(
        "astrocyte_feature_extraction.segmentation._cellpose_version",
        lambda: "test",
    )
    model = FakeModel()
    mask_path, metadata = segment_record(record, config, tmp_path / "output", model=model)
    labels = tifffile.imread(mask_path)
    assert set(np.unique(labels)) == {0, 1, 2}
    assert metadata["cell_count"] == 2
    assert metadata["settings"]["input_channel_indices"] == [1, 2, 0]
    assert model.input_image is not None
    assert np.all(model.input_image[..., 0] == 1000)
    assert np.all(model.input_image[..., 1] == 2000)
    assert np.array_equal(
        model.input_image[..., 2],
        np.arange(record.height * record.width, dtype=np.uint16).reshape(
            record.height, record.width
        ),
    )

    table = export_cell_crops(
        record,
        mask_path,
        tmp_path / "output",
        {
            "margin_pixels": 2,
            "write_article_64x64": True,
            "article_size": 64,
            "article_epsilon": 1.0e-5,
        },
    )
    assert list(table["cell_id"]) == [1, 2]
    crop = np.moveaxis(tifffile.imread(table.iloc[0]["crop_path"]), 0, -1)
    crop_mask = tifffile.imread(table.iloc[0]["crop_mask_path"]).astype(bool)
    assert np.all(crop[~crop_mask] == 0)
    normalized = np.moveaxis(
        tifffile.imread(table.iloc[0]["article_crop_path"]), 0, -1
    )
    assert normalized.shape == (64, 64, 3)
    assert np.all(np.isfinite(normalized))


def test_article_preprocess_is_channelwise() -> None:
    crop = np.zeros((8, 12, 3), dtype=np.uint8)
    crop[2:6, 3:9, 0] = 100
    crop[2:6, 3:9, 1] = 200
    normalized = article_preprocess(crop)
    assert normalized.shape == (64, 64, 3)
    assert abs(float(normalized[..., 0].mean())) < 1.0e-5
    assert abs(float(normalized[..., 1].mean())) < 1.0e-5
    assert np.all(normalized[..., 2] == 0)
