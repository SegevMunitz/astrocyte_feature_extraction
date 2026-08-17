"""Microscopy image loading helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from .inventory import ImageRecord


def load_hwc(record: ImageRecord) -> np.ndarray:
    """Load one static multichannel image in height-width-channel order."""
    array = np.asarray(tifffile.imread(record.path))
    if tuple(array.shape) != tuple(record.shape):
        raise ValueError(
            f"Image shape changed after inventory: {record.path}: "
            f"expected {record.shape}, got {array.shape}"
        )
    array = np.moveaxis(array, record.channel_axis, -1)
    if array.ndim != 3 or array.shape[-1] != record.channel_count:
        raise ValueError(f"Expected HWC array after channel move, got {array.shape}")
    return array


def write_hwc(path: Path, array: np.ndarray) -> None:
    """Write a multichannel TIFF with explicit YXC axes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        np.moveaxis(array, -1, 0),
        metadata={"axes": "CYX"},
        photometric="minisblack",
    )
