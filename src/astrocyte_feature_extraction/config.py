"""Configuration loading and validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when pipeline configuration is incomplete or inconsistent."""


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration and resolve environment variables."""
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
    return config


def output_path(config: dict[str, Any], *parts: str) -> Path:
    """Return a path below the configured output root."""
    return Path(config["paths"]["output_dir"]).expanduser().joinpath(*parts)
