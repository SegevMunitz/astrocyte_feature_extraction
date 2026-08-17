"""Cellpose inference with explicit three-channel settings and provenance."""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import json
from dataclasses import asdict
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from skimage.segmentation import find_boundaries

from .images import load_hwc
from .inventory import ImageRecord


class SegmentationError(RuntimeError):
    """Raised when Cellpose cannot safely produce an instance mask."""


@lru_cache(maxsize=8)
def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute a streaming SHA-256 checksum."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _cellpose_version() -> str:
    try:
        return importlib.metadata.version("cellpose")
    except importlib.metadata.PackageNotFoundError as exc:
        raise SegmentationError(
            "Cellpose is not installed. Install the project with the segmentation extra."
        ) from exc


def _build_model(checkpoint: Path, gpu: bool, channel_count: int) -> Any:
    try:
        from cellpose import models
    except ImportError as exc:
        raise SegmentationError(
            "Cellpose is not installed. Install the project with the segmentation extra."
        ) from exc
    kwargs: dict[str, Any] = {"gpu": gpu, "pretrained_model": str(checkpoint)}
    signature = inspect.signature(models.CellposeModel)
    if "nchan" in signature.parameters:
        kwargs["nchan"] = channel_count
    try:
        return models.CellposeModel(**kwargs)
    except Exception as exc:  # Cellpose wraps backend/device errors inconsistently.
        raise SegmentationError(f"Could not load Cellpose checkpoint {checkpoint}: {exc}") from exc


def _to_uint8_plane(plane: np.ndarray) -> np.ndarray:
    finite = np.asarray(plane, dtype=np.float32)
    low, high = np.nanpercentile(finite, (1.0, 99.0))
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return np.zeros(plane.shape, dtype=np.uint8)
    return np.clip((finite - low) * (255.0 / (high - low)), 0, 255).astype(np.uint8)


def _write_overlay(path: Path, image: np.ndarray, labels: np.ndarray) -> None:
    base = _to_uint8_plane(np.max(image, axis=-1))
    overlay = np.repeat(base[..., None], 3, axis=-1)
    overlay[find_boundaries(labels, mode="outer")] = np.array([255, 0, 0], dtype=np.uint8)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, overlay, photometric="rgb")


def _eval_kwargs(model: Any, settings: dict[str, Any]) -> dict[str, Any]:
    requested: dict[str, Any] = {
        "channels": None,
        "channel_axis": -1,
        "diameter": settings.get("diameter"),
        "flow_threshold": float(settings["flow_threshold"]),
        "cellprob_threshold": float(settings["cellprob_threshold"]),
        "normalize": bool(settings["normalize"]),
        "min_size": int(settings["min_size"]),
    }
    signature = inspect.signature(model.eval)
    return {key: value for key, value in requested.items() if key in signature.parameters}


def _input_channel_indices(settings: dict[str, Any], channel_count: int) -> tuple[int, ...]:
    requested = settings.get("input_channel_indices")
    if requested is None:
        return tuple(range(channel_count))
    if not isinstance(requested, (list, tuple)) or not requested:
        raise SegmentationError("cellpose.input_channel_indices must be a non-empty list")
    indices = tuple(int(index) for index in requested)
    if any(index < 0 or index >= channel_count for index in indices):
        raise SegmentationError(
            "cellpose.input_channel_indices contains an index outside the source image channels"
        )
    return indices


def _provenance_settings(
    model: Any, settings: dict[str, Any], channel_count: int
) -> dict[str, Any]:
    return {
        "eval": _eval_kwargs(model, settings),
        "input_channel_indices": list(_input_channel_indices(settings, channel_count)),
    }


def segment_record(
    record: ImageRecord,
    config: dict[str, Any],
    output_root: Path,
    model: Any | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Segment one image and write a lossless integer label image."""
    checkpoint = Path(config["paths"]["checkpoint"]).expanduser()
    if not checkpoint.is_file():
        raise SegmentationError(f"Cellpose checkpoint does not exist: {checkpoint}")
    settings = config["cellpose"]
    installed_version = _cellpose_version()
    expected_version = settings.get("expected_version")
    if expected_version and installed_version != str(expected_version):
        raise SegmentationError(
            f"Cellpose version mismatch: expected {expected_version}, found {installed_version}"
        )
    image = load_hwc(record)
    input_channel_indices = _input_channel_indices(settings, record.channel_count)
    model_image = image[..., input_channel_indices]
    if model is None:
        model = _build_model(checkpoint, bool(settings["gpu"]), model_image.shape[-1])
    eval_kwargs = _eval_kwargs(model, settings)
    try:
        result = model.eval(model_image, **eval_kwargs)
    except Exception as exc:
        raise SegmentationError(f"Cellpose inference failed for {record.path}: {exc}") from exc
    masks = result[0] if isinstance(result, tuple) else result
    labels = np.asarray(masks)
    if labels.shape != (record.height, record.width):
        raise SegmentationError(
            f"Cellpose returned shape {labels.shape}; expected {(record.height, record.width)}"
        )
    if not np.issubdtype(labels.dtype, np.integer):
        if not np.all(np.equal(labels, np.floor(labels))):
            raise SegmentationError("Cellpose returned non-integral labels")
        labels = labels.astype(np.uint32)
    if np.any(labels < 0):
        raise SegmentationError("Cellpose returned negative labels")
    labels = labels.astype(np.uint32, copy=False)
    unique = np.unique(labels)
    positive = unique[unique > 0]
    if positive.size and not np.array_equal(positive, np.arange(1, positive.size + 1)):
        remapped = np.zeros_like(labels, dtype=np.uint32)
        for new_id, old_id in enumerate(positive, start=1):
            remapped[labels == old_id] = new_id
        labels = remapped
    mask_path = output_root / "masks" / f"{record.image_id}_labels.tif"
    mask_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(mask_path, labels, metadata={"axes": "YX"}, photometric="minisblack")
    _write_overlay(output_root / "qc" / f"{record.image_id}_overlay.tif", image, labels)
    border_ids = np.unique(
        np.concatenate((labels[0], labels[-1], labels[:, 0], labels[:, -1]))
    )
    border_ids = border_ids[border_ids > 0]
    metadata = {
        "image": asdict(record),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256_file(checkpoint),
        "cellpose_version": installed_version,
        "settings": _provenance_settings(model, settings, record.channel_count),
        "mask_path": str(mask_path),
        "cell_count": int(labels.max()),
        "empty_segmentation": bool(labels.max() == 0),
        "border_cell_ids": [int(value) for value in border_ids],
    }
    manifest_path = output_root / "manifests" / f"{record.image_id}_segmentation.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    return mask_path, metadata


def segment_all(
    records: list[ImageRecord], config: dict[str, Any], output_root: Path
) -> dict[str, Path]:
    """Load the model once and segment every inventory record."""
    if not records:
        return {}
    checkpoint = Path(config["paths"]["checkpoint"]).expanduser()
    model = _build_model(
        checkpoint,
        bool(config["cellpose"]["gpu"]),
        len(_input_channel_indices(config["cellpose"], records[0].channel_count)),
    )
    expected_settings = _provenance_settings(
        model, config["cellpose"], records[0].channel_count
    )
    expected_checkpoint_hash = sha256_file(checkpoint)
    expected_version = _cellpose_version()
    pinned_version = config["cellpose"].get("expected_version")
    if pinned_version and expected_version != str(pinned_version):
        raise SegmentationError(
            f"Cellpose version mismatch: expected {pinned_version}, found {expected_version}"
        )
    outputs: dict[str, Path] = {}
    for record in records:
        if record.channel_count != records[0].channel_count:
            raise SegmentationError("All images must have the same channel count")
        mask_path = output_root / "masks" / f"{record.image_id}_labels.tif"
        manifest_path = output_root / "manifests" / f"{record.image_id}_segmentation.json"
        if (
            config["cellpose"].get("resume", True)
            and mask_path.is_file()
            and manifest_path.is_file()
        ):
            try:
                metadata = json.loads(manifest_path.read_text(encoding="utf-8"))
                mask_shape = tuple(tifffile.imread(mask_path).shape)
                is_compatible = (
                    metadata.get("checkpoint_sha256") == expected_checkpoint_hash
                    and metadata.get("cellpose_version") == expected_version
                    and metadata.get("settings") == expected_settings
                    and mask_shape == (record.height, record.width)
                )
            except (OSError, ValueError, json.JSONDecodeError):
                is_compatible = False
            if is_compatible:
                outputs[record.image_id] = mask_path
                continue
        mask_path, _ = segment_record(record, config, output_root, model=model)
        outputs[record.image_id] = mask_path
    return outputs
