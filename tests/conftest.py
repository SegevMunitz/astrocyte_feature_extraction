from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from astrocyte_feature_extraction.inventory import ImageRecord


def write_test_image(path: Path, height: int = 12, width: int = 16) -> ImageRecord:
    channels = np.zeros((3, height, width), dtype=np.uint16)
    channels[0] = np.arange(height * width, dtype=np.uint16).reshape(height, width)
    channels[1] = 1000
    channels[2] = 2000
    tifffile.imwrite(
        path,
        channels,
        metadata={"axes": "CYX"},
        photometric="minisblack",
    )
    return ImageRecord(
        image_id="image_a",
        path=str(path),
        relative_path=path.name,
        shape=channels.shape,
        dtype=str(channels.dtype),
        axes="CYX",
        channel_axis=0,
        channel_count=3,
        height=height,
        width=width,
        pixel_size_x_um=0.5,
        pixel_size_y_um=0.5,
    )
