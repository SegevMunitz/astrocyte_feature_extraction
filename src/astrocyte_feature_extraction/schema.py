"""CellProfiler CSV schema and pipeline provenance inspection."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


class SchemaError(RuntimeError):
    """Raised when the requested measurement schema is ambiguous."""


CATEGORY_MODULES: dict[str, str] = {
    "AreaShape": "MeasureObjectSizeShape",
    "Intensity": "MeasureObjectIntensity",
    "Texture": "MeasureTexture",
    "Neighbors": "MeasureObjectNeighbors",
    "RadialDistribution": "MeasureObjectRadialDistribution",
    "Granularity": "MeasureGranularity",
    "Correlation": "MeasureColocalization",
    "ObjectSkeleton": "MeasureObjectSkeleton",
    "Location": "IdentifyObjects/MeasureObjectIntensity",
    "Number": "IdentifyObjects",
    "Parent": "RelateObjects",
    "Children": "RelateObjects",
}

IDENTIFIER_PREFIXES = (
    "ImageNumber",
    "ObjectNumber",
    "Metadata_",
    "FileName_",
    "PathName_",
    "URL_",
    "Group_",
    "ExecutionTime_",
    "ModuleError_",
)


@dataclass(frozen=True)
class MetricSpec:
    column: str
    category: str
    module: str
    role: str


@dataclass(frozen=True)
class PipelineSpec:
    path: str
    pipeline_format_version: str | None
    date_revision: str | None
    git_hash: str | None
    modules: tuple[str, ...]


def read_csv_header(path: str | Path) -> list[str]:
    """Read the first CSV record using pandas-compatible duplicate names."""
    csv_path = Path(path).expanduser()
    if not csv_path.is_file():
        raise SchemaError(f"CellProfiler CSV does not exist: {csv_path}")
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        try:
            header = next(csv.reader(handle))
        except StopIteration as exc:
            raise SchemaError(f"CellProfiler CSV is empty: {csv_path}") from exc
    cleaned = [column.strip() for column in header]
    while cleaned and not cleaned[-1]:
        cleaned.pop()
    if not cleaned:
        raise SchemaError("CSV header contains no columns")
    if any(not column for column in cleaned):
        raise SchemaError("CSV header contains an empty non-trailing column name")
    occurrence_counts: dict[str, int] = {}
    canonical: list[str] = []
    for column in cleaned:
        occurrence = occurrence_counts.get(column, 0)
        canonical.append(column if occurrence == 0 else f"{column}.{occurrence}")
        occurrence_counts[column] = occurrence + 1
    return canonical


def classify_column(column: str) -> MetricSpec:
    """Map a CellProfiler export column to its producing module family."""
    if column.startswith(IDENTIFIER_PREFIXES):
        return MetricSpec(column, "Identifier", "CellProfiler export/input", "provenance")
    tokens = column.split("_")
    for category, module in CATEGORY_MODULES.items():
        if category in tokens or column.startswith(f"{category}_"):
            return MetricSpec(column, category, module, "feature")
    if column in {"ImageNumber", "ObjectNumber"}:
        return MetricSpec(column, "Identifier", "CellProfiler export/input", "provenance")
    return MetricSpec(column, "Unmapped", "Unmapped", "feature")


def build_metric_manifest(columns: list[str], fail_on_unmapped: bool = True) -> list[MetricSpec]:
    """Create an ordered metric manifest from exported CSV columns."""
    specs = [classify_column(column) for column in columns]
    unmapped = [spec.column for spec in specs if spec.category == "Unmapped"]
    if fail_on_unmapped and unmapped:
        raise SchemaError(
            "Unmapped CellProfiler columns require an explicit module mapping: "
            + ", ".join(unmapped)
        )
    return specs


def inspect_cppipe(path: Path) -> PipelineSpec:
    """Inspect a text CellProfiler pipeline without executing it."""
    text = path.read_text(encoding="utf-8-sig", errors="replace")

    def header_value(name: str) -> str | None:
        match = re.search(rf"^{re.escape(name)}:(.*)$", text, flags=re.MULTILINE)
        return match.group(1).strip() if match else None

    modules = tuple(
        match.group(1)
        for match in re.finditer(
            r"^([A-Za-z][A-Za-z0-9]+):\[module_num:\d+\|",
            text,
            flags=re.MULTILINE,
        )
    )
    return PipelineSpec(
        path=str(path),
        pipeline_format_version=header_value("Version"),
        date_revision=header_value("DateRevision"),
        git_hash=header_value("GitHash"),
        modules=modules,
    )


def discover_pipelines(directory: str | Path) -> list[PipelineSpec]:
    """Find and inspect all text pipelines under a directory."""
    root = Path(directory).expanduser()
    if not root.is_dir():
        raise SchemaError(f"CellProfiler pipeline directory does not exist: {root}")
    pipelines = [inspect_cppipe(path) for path in sorted(root.rglob("*.cppipe"))]
    if not pipelines:
        projects = list(root.rglob("*.cpproj"))
        detail = f" Found {len(projects)} .cpproj files but no inspectable .cppipe." if projects else ""
        raise SchemaError(f"No .cppipe files found below {root}.{detail}")
    return pipelines


def validate_pipeline_modules(
    metrics: list[MetricSpec], pipelines: list[PipelineSpec]
) -> dict[str, list[str]]:
    """Return required measurement modules and candidate pipelines containing each."""
    required = sorted(
        {
            spec.module
            for spec in metrics
            if spec.role == "feature" and spec.module not in {"Unmapped", "IdentifyObjects/MeasureObjectIntensity"}
        }
    )
    return {
        module: [pipeline.path for pipeline in pipelines if module in pipeline.modules]
        for module in required
    }


def write_schema_manifest(
    metrics: list[MetricSpec],
    pipelines: list[PipelineSpec],
    destination: Path,
) -> None:
    """Write requested measurements and pipeline provenance to JSON."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "metrics": [asdict(metric) for metric in metrics],
        "pipelines": [
            {**asdict(pipeline), "modules": list(pipeline.modules)} for pipeline in pipelines
        ],
        "module_candidates": validate_pipeline_modules(metrics, pipelines),
    }
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
