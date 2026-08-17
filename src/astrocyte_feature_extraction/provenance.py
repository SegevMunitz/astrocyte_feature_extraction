"""Run-level reproducibility manifests."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
from typing import Any

from .inventory import ImageRecord


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def write_run_manifest(
    config: dict[str, Any],
    records: list[ImageRecord],
    output_root: Path,
    config_path: Path,
) -> Path:
    """Record configuration, input hashes, and relevant package versions."""
    inputs = []
    for record in records:
        path = Path(record.path)
        stat = path.stat()
        inputs.append(
            {
                "image_id": record.image_id,
                "path": str(path),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": _sha256(path),
            }
        )
    payload = {
        "config_path": str(config_path.resolve()),
        "config": config,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": {
            name: _version(name)
            for name in (
                "astrocyte-feature-extraction",
                "cellpose",
                "cellprofiler",
                "numpy",
                "pandas",
                "scipy",
                "scikit-image",
                "tifffile",
            )
        },
        "inputs": inputs,
    }
    destination = output_root / "manifests" / "run_manifest.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    return destination
