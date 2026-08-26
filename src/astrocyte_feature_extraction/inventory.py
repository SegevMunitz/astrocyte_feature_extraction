"""Discover microscopy inputs and capture their dimensional metadata."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import tifffile


class InventoryError(RuntimeError):
    """Raised when image metadata cannot be interpreted safely."""


@dataclass(frozen=True)
class ImageRecord:
    image_id: str
    path: str
    relative_path: str
    shape: tuple[int, ...]
    dtype: str
    axes: str | None
    channel_axis: int
    channel_count: int
    height: int
    width: int
    pixel_size_x_um: float | None
    pixel_size_y_um: float | None


def _stable_image_id(relative_path: Path) -> str:
    readable = re.sub(r"[^A-Za-z0-9._-]+", "_", relative_path.as_posix())
    digest = hashlib.sha1(relative_path.as_posix().encode("utf-8")).hexdigest()[:10]
    return f"{Path(readable).stem}_{digest}"


def _ome_pixel_size(ome_xml: str | None) -> tuple[float | None, float | None]:
    if not ome_xml:
        return None, None
    values: list[float | None] = []
    for key in ("PhysicalSizeX", "PhysicalSizeY"):
        match = re.search(rf'{key}="([^"]+)"', ome_xml)
        values.append(float(match.group(1)) if match else None)
    return values[0], values[1]


def _resolve_channel_axis(
    shape: tuple[int, ...],
    axes: str | None,
    configured_axis: str | int,
    expected_channels: int,
) -> int:
    if configured_axis != "auto":
        axis = int(configured_axis)
        if axis < 0:
            axis += len(shape)
        if axis < 0 or axis >= len(shape):
            raise InventoryError(f"Channel axis {configured_axis} is invalid for shape {shape}")
        if shape[axis] != expected_channels:
            raise InventoryError(
                f"Configured channel axis has {shape[axis]} channels; expected {expected_channels}"
            )
        return axis
    if axes and "C" in axes:
        axis = axes.index("C")
        if shape[axis] != expected_channels:
            raise InventoryError(
                f"OME channel axis has {shape[axis]} channels; expected {expected_channels}"
            )
        return axis
    candidates = [index for index, size in enumerate(shape) if size == expected_channels]
    if len(candidates) != 1:
        raise InventoryError(
            f"Cannot infer one channel axis for shape {shape}; candidates={candidates}. "
            "Set images.channel_axis explicitly."
        )
    return candidates[0]


def inspect_image(
    path: Path,
    root: Path,
    configured_axis: str | int,
    expected_channels: int,
    configured_pixel_size_um: float | None = None,
) -> ImageRecord:
    """Read TIFF metadata without loading the full pixel array."""
    try:
        with tifffile.TiffFile(path) as tif:
            series = tif.series[0]
            shape = tuple(int(value) for value in series.shape)
            axes = getattr(series, "axes", None)
            dtype = str(series.dtype)
            px_x, px_y = _ome_pixel_size(tif.ome_metadata)
    except (OSError, ValueError, tifffile.TiffFileError) as exc:
        raise InventoryError(f"Could not read TIFF metadata for {path}: {exc}") from exc
    channel_axis = _resolve_channel_axis(shape, axes, configured_axis, expected_channels)
    spatial_axes = [index for index in range(len(shape)) if index != channel_axis]
    if len(spatial_axes) != 2:
        raise InventoryError(
            f"Only static 2D multichannel images are supported; got shape={shape}, axes={axes}"
        )
    if configured_pixel_size_um is not None:
        px_x = px_y = float(configured_pixel_size_um)
    relative = path.relative_to(root)
    return ImageRecord(
        image_id=_stable_image_id(relative),
        path=str(path),
        relative_path=relative.as_posix(),
        shape=shape,
        dtype=dtype,
        axes=axes,
        channel_axis=channel_axis,
        channel_count=shape[channel_axis],
        height=shape[spatial_axes[-2]],
        width=shape[spatial_axes[-1]],
        pixel_size_x_um=px_x,
        pixel_size_y_um=px_y,
    )


def discover_images(config: dict[str, Any]) -> list[ImageRecord]:
    """Discover and validate all configured static multichannel TIFF images."""
    image_config = config["images"]
    root = Path(config["paths"]["images_dir"]).expanduser()
    if not root.is_dir():
        raise InventoryError(f"Image directory does not exist: {root}")
    suffixes = {suffix.lower() for suffix in image_config["extensions"]}
    include_pattern = image_config.get("include_regex")
    include_regex = re.compile(str(include_pattern)) if include_pattern else None
    iterator: Iterable[Path] = root.rglob("*") if image_config.get("recursive", True) else root.glob("*")
    output_root = Path(config["paths"]["output_dir"]).expanduser()
    candidates = sorted(
        path
        for path in iterator
        if path.is_file()
        and any(path.name.lower().endswith(suffix) for suffix in suffixes)
        and output_root not in path.parents
        and "cellprofiler results" not in {part.lower() for part in path.parts}
        and (include_regex is None or include_regex.search(path.name))
    )
    maximum_images = image_config.get("max_images")
    if maximum_images is not None:
        candidates = candidates[: int(maximum_images)]
    if not candidates:
        raise InventoryError(f"No configured TIFF images found below {root}")
    skip_incompatible = bool(image_config.get("skip_incompatible", False))
    records: list[ImageRecord] = []
    skipped: list[dict[str, str]] = []
    for path in candidates:
        try:
            records.append(
                inspect_image(
                    path,
                    root,
                    image_config.get("channel_axis", "auto"),
                    int(image_config["expected_channels"]),
                    image_config.get("pixel_size_um"),
                )
            )
        except InventoryError as exc:
            if not skip_incompatible:
                raise
            skipped.append({"path": str(path), "reason": str(exc)})
    if skip_incompatible and skipped:
        skip_path = output_root / "manifests" / "skipped_images.json"
        skip_path.parent.mkdir(parents=True, exist_ok=True)
        skip_path.write_text(json.dumps(skipped, indent=2, sort_keys=True), encoding="utf-8")
    if not records:
        raise InventoryError(
            f"No compatible TIFF images found below {root}"
            + (f" ({len(skipped)} skipped)" if skipped else "")
        )
    return records


def write_inventory(records: list[ImageRecord], destination: Path) -> None:
    """Write a deterministic JSON inventory."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps([asdict(record) for record in records], indent=2, sort_keys=True),
        encoding="utf-8",
    )
