from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from astrocyte_feature_extraction.spatial import knn_edges, nodes_from_labels

from conftest import write_test_image


def test_nodes_and_knn_are_global_and_same_image(tmp_path: Path) -> None:
    record = write_test_image(tmp_path / "image.tif", height=20, width=30)
    labels = np.zeros((20, 30), dtype=np.uint16)
    labels[1:4, 1:4] = 1
    labels[8:11, 10:13] = 2
    labels[15:18, 25:28] = 3
    labels_path = tmp_path / "labels.tif"
    tifffile.imwrite(labels_path, labels, metadata={"axes": "YX"})

    nodes = nodes_from_labels(record, labels_path)
    assert list(nodes["cell_id"]) == [1, 2, 3]
    assert nodes["centroid_x_normalized"].between(0, 1).all()
    assert nodes.loc[nodes["cell_id"] == 1, "centroid_x_um"].item() == 1.0

    edges = knn_edges(nodes, k=1)
    assert len(edges) == 3
    assert set(edges["image_id"]) == {"image_a"}
    assert (edges["source_cell_id"] != edges["target_cell_id"]).all()
    assert (edges["neighbor_rank"] == 1).all()


def test_node_qc_flags_disconnected_label(tmp_path: Path) -> None:
    record = write_test_image(tmp_path / "image.tif", height=20, width=30)
    labels = np.zeros((20, 30), dtype=np.uint16)
    labels[2:4, 2:4] = 1
    labels[10:12, 10:12] = 1
    labels_path = tmp_path / "labels.tif"
    tifffile.imwrite(labels_path, labels, metadata={"axes": "YX"})

    nodes = nodes_from_labels(
        record,
        labels_path,
        {"min_cell_area_pixels": 10, "max_cell_area_pixels": 100},
    )
    assert nodes.iloc[0]["fragment_count"] == 2
    assert bool(nodes.iloc[0]["is_fragmented"])
    assert bool(nodes.iloc[0]["is_too_small"])
    assert not bool(nodes.iloc[0]["is_too_large"])
