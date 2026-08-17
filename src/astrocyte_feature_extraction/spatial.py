"""Per-cell global coordinates and graph-ready neighborhood tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from scipy.spatial import cKDTree

from .inventory import ImageRecord


class SpatialError(RuntimeError):
    """Raised when label geometry cannot produce valid spatial features."""


def nodes_from_labels(
    record: ImageRecord,
    labels_path: Path,
    qc_config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Compute global object centroids and border quality fields."""
    labels = np.asarray(tifffile.imread(labels_path))
    if labels.shape != (record.height, record.width):
        raise SpatialError(
            f"Label dimensions changed for {record.image_id}: {labels.shape}"
        )
    cell_ids = np.unique(labels)
    cell_ids = cell_ids[cell_ids > 0].astype(int)
    rows: list[dict[str, Any]] = []
    qc_config = qc_config or {}
    minimum_area = int(qc_config.get("min_cell_area_pixels", 0))
    maximum_value = qc_config.get("max_cell_area_pixels")
    maximum_area = int(maximum_value) if maximum_value is not None else None
    for cell_id in cell_ids:
        selected = labels == cell_id
        y, x = ndimage.center_of_mass(selected)
        yy, xx = np.nonzero(selected)
        fragment_count = int(
            ndimage.label(selected, structure=np.ones((3, 3), dtype=np.uint8))[1]
        )
        border_distance = min(
            float(x),
            float(y),
            float(record.width - 1 - x),
            float(record.height - 1 - y),
        )
        touches_border = bool(
            np.any(yy == 0)
            or np.any(xx == 0)
            or np.any(yy == record.height - 1)
            or np.any(xx == record.width - 1)
        )
        row: dict[str, Any] = {
            "image_id": record.image_id,
            "cell_id": int(cell_id),
            "global_cell_id": f"{record.image_id}:{int(cell_id)}",
            "centroid_x_px": float(x),
            "centroid_y_px": float(y),
            "centroid_x_normalized": float(x / max(record.width - 1, 1)),
            "centroid_y_normalized": float(y / max(record.height - 1, 1)),
            "distance_to_border_px": border_distance,
            "touches_image_border": touches_border,
            "area_pixels": int(yy.size),
            "fragment_count": fragment_count,
            "is_fragmented": fragment_count > 1,
            "is_too_small": int(yy.size) < minimum_area,
            "is_too_large": maximum_area is not None and int(yy.size) > maximum_area,
        }
        if record.pixel_size_x_um is not None and record.pixel_size_y_um is not None:
            row.update(
                {
                    "centroid_x_um": float(x * record.pixel_size_x_um),
                    "centroid_y_um": float(y * record.pixel_size_y_um),
                    "distance_to_border_um": float(
                        min(
                            x * record.pixel_size_x_um,
                            y * record.pixel_size_y_um,
                            (record.width - 1 - x) * record.pixel_size_x_um,
                            (record.height - 1 - y) * record.pixel_size_y_um,
                        )
                    ),
                }
            )
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def knn_edges(nodes: pd.DataFrame, k: int) -> pd.DataFrame:
    """Create directed same-image k-nearest-neighbor edges."""
    if k < 1:
        raise SpatialError("spatial.k_neighbors must be at least 1")
    rows: list[dict[str, Any]] = []
    for image_id, group in nodes.groupby("image_id", sort=True):
        ordered = group.sort_values("cell_id").reset_index(drop=True)
        count = len(ordered)
        if count < 2:
            continue
        coordinates = ordered[["centroid_x_px", "centroid_y_px"]].to_numpy(float)
        has_physical_coordinates = (
            {"centroid_x_um", "centroid_y_um"}.issubset(ordered.columns)
            and ordered[["centroid_x_um", "centroid_y_um"]].notna().all().all()
        )
        neighbor_count = min(k, count - 1)
        distances, indices = cKDTree(coordinates).query(
            coordinates, k=neighbor_count + 1
        )
        if distances.ndim == 1:
            distances = distances[:, None]
            indices = indices[:, None]
        for source_index in range(count):
            source = ordered.iloc[source_index]
            candidates = [
                (float(distance), int(target_index))
                for distance, target_index in zip(
                    distances[source_index], indices[source_index], strict=True
                )
                if int(target_index) != source_index
            ]
            candidates.sort(key=lambda item: (item[0], int(ordered.iloc[item[1]]["cell_id"])))
            for rank, (distance_px, target_index) in enumerate(
                candidates[:neighbor_count], start=1
            ):
                target = ordered.iloc[target_index]
                dx_px = float(target["centroid_x_px"] - source["centroid_x_px"])
                dy_px = float(target["centroid_y_px"] - source["centroid_y_px"])
                edge: dict[str, Any] = {
                    "image_id": image_id,
                    "source_cell_id": int(source["cell_id"]),
                    "target_cell_id": int(target["cell_id"]),
                    "source_global_cell_id": source["global_cell_id"],
                    "target_global_cell_id": target["global_cell_id"],
                    "neighbor_rank": rank,
                    "distance_px": distance_px,
                    "delta_x_px": dx_px,
                    "delta_y_px": dy_px,
                }
                if has_physical_coordinates:
                    dx_um = float(target["centroid_x_um"] - source["centroid_x_um"])
                    dy_um = float(target["centroid_y_um"] - source["centroid_y_um"])
                    edge.update(
                        {
                            "distance_um": float(np.hypot(dx_um, dy_um)),
                            "delta_x_um": dx_um,
                            "delta_y_um": dy_um,
                        }
                    )
                rows.append(edge)
    columns = [
        "image_id",
        "source_cell_id",
        "target_cell_id",
        "source_global_cell_id",
        "target_global_cell_id",
        "neighbor_rank",
        "distance_px",
        "delta_x_px",
        "delta_y_px",
    ]
    optional_physical = ["distance_um", "delta_x_um", "delta_y_um"]
    if any("distance_um" in row for row in rows):
        columns.extend(optional_physical)
    return pd.DataFrame.from_records(rows, columns=columns)


def build_spatial_tables(
    records: list[ImageRecord],
    mask_paths: dict[str, Path],
    output_root: Path,
    spatial_config: dict[str, Any],
    qc_config: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build and save node and edge tables for all images."""
    node_tables = [
        nodes_from_labels(record, mask_paths[record.image_id], qc_config)
        for record in records
    ]
    nodes = pd.concat(node_tables, ignore_index=True) if node_tables else pd.DataFrame()
    edges = knn_edges(nodes, int(spatial_config["k_neighbors"])) if not nodes.empty else pd.DataFrame()
    table_dir = output_root / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(table_dir / "nodes.csv", index=False)
    nodes.to_parquet(table_dir / "nodes.parquet", index=False)
    edges.to_csv(table_dir / "edges.csv", index=False)
    edges.to_parquet(table_dir / "edges.parquet", index=False)
    return nodes, edges
