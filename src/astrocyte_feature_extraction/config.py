"""Configuration loading and validation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when pipeline configuration is incomplete or inconsistent."""


_PATH_KEYS = (
    "images_dir",
    "checkpoint",
    "cellprofiler_results_csv",
    "cellprofiler_pipeline_dir",
    "cellprofiler_measurement_pipeline",
    "output_dir",
    "bundle_dir",
    "cluster_assigner",
    "cluster_csv",
    "morph_csv",
    "crop_measure_script",
    "slurm_log_dir",
)


def _resolve_path(value: str, root: Path) -> str:
    path = Path(value).expanduser()
    if path.is_absolute():
        return str(path)
    return str((root / path).resolve())


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration and resolve paths relative to ROOT or config dir."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        raise ConfigError(f"Configuration file does not exist: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ConfigError("Configuration root must be a mapping")
    for section in (
        "paths",
        "images",
        "cellpose",
        "crops",
        "spatial",
        "qc",
        "cellprofiler",
    ):
        if section not in config or not isinstance(config[section], dict):
            raise ConfigError(f"Missing configuration section: {section}")
    required_paths = (
        "images_dir",
        "checkpoint",
        "cellprofiler_results_csv",
        "cellprofiler_pipeline_dir",
        "output_dir",
    )
    for key in required_paths:
        if not config["paths"].get(key):
            raise ConfigError(f"Missing paths.{key}")

    root = Path(os.environ.get("ROOT", "")).expanduser() if os.environ.get("ROOT") else None
    if root is None or not root.is_dir():
        # Prefer explicit paths.root, else directory containing configs/
        explicit = config["paths"].get("root")
        if explicit:
            root = Path(explicit).expanduser().resolve()
        else:
            root = config_path.parent.parent if config_path.parent.name == "configs" else config_path.parent

    for key in _PATH_KEYS:
        value = config["paths"].get(key)
        if value:
            config["paths"][key] = _resolve_path(str(value), root)

    for key in ("executable", "python_executable", "direct_script"):
        value = config["cellprofiler"].get(key)
        if value and not Path(str(value)).expanduser().is_absolute():
            config["cellprofiler"][key] = _resolve_path(str(value), root)

    config["paths"]["root"] = str(root)
    return config


def output_path(config: dict[str, Any], *parts: str) -> Path:
    """Return a path below the configured output root."""
    return Path(config["paths"]["output_dir"]).expanduser().joinpath(*parts)
