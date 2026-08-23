"""Per-cell PCA and plots on a full single-cell feature table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


META_PREFIXES = (
    "sample",
    "crop_path",
    "mask_path",
    "Time",
    "PhotoID",
    "cluster",
    "original_cell_id",
    "compact_cell_id",
    "measurement_backend",
    "crop_size_pixels",
)


def _feature_columns(table: pd.DataFrame) -> list[str]:
    columns: list[str] = []
    for column in table.columns:
        if column in META_PREFIXES or column.startswith("Metadata_"):
            continue
        if pd.api.types.is_numeric_dtype(table[column]):
            columns.append(column)
    return columns


def _parse_time(sample: pd.Series) -> pd.Series:
    return sample.astype(str).str.extract(r"^(ctrl|4h|24h|72h|7d)", expand=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/"
            "results/single_cell_full_features.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/"
            "results/plots/per_cell_pca"
        ),
    )
    parser.add_argument("--n-components", type=int, default=10)
    args = parser.parse_args()

    table = pd.read_csv(args.input)
    if "Time" not in table.columns and "sample" in table.columns:
        table["Time"] = _parse_time(table["sample"])

    features = _feature_columns(table)
    matrix = table[features].replace([np.inf, -np.inf], np.nan)
    keep = matrix.columns[matrix.notna().any(axis=0)]
    matrix = matrix[keep].fillna(matrix[keep].median(numeric_only=True))
    variance = matrix.var(axis=0, numeric_only=True)
    keep = variance[variance > 1.0e-12].index.tolist()
    matrix = matrix[keep]

    scaled = StandardScaler().fit_transform(matrix)
    n_components = min(args.n_components, scaled.shape[1], scaled.shape[0])
    pca = PCA(n_components=n_components, random_state=0)
    scores = pca.fit_transform(scaled)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    score_table = table[
        [column for column in ("sample", "Time", "cluster", "original_cell_id") if column in table.columns]
    ].copy()
    for index in range(n_components):
        score_table[f"PC{index + 1}"] = scores[:, index]
    score_table.to_csv(output_dir / "pca_scores.csv", index=False)

    loadings = pd.DataFrame(
        pca.components_.T,
        index=keep,
        columns=[f"PC{index + 1}" for index in range(n_components)],
    )
    loadings.to_csv(output_dir / "pca_loadings.csv")

    explained = pd.DataFrame(
        {
            "component": [f"PC{index + 1}" for index in range(n_components)],
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    explained.to_csv(output_dir / "pca_explained_variance.csv", index=False)

    time_colors = {
        "ctrl": "#999999",
        "4h": "#E41A1C",
        "24h": "#377EB8",
        "72h": "#4DAF4A",
        "7d": "#984EA3",
    }

    fig, axis = plt.subplots(figsize=(8, 6))
    if "Time" in table.columns and table["Time"].notna().any():
        for time_point, group in score_table.groupby("Time", dropna=True):
            axis.scatter(
                group["PC1"],
                group["PC2"],
                s=8,
                alpha=0.45,
                label=str(time_point),
                c=time_colors.get(str(time_point), "#555555"),
            )
        axis.legend(title="Time", frameon=False)
    else:
        axis.scatter(scores[:, 0], scores[:, 1], s=8, alpha=0.45, c="#555555")
    axis.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    axis.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    axis.set_title("Per-cell PCA on full CellProfiler features")
    fig.tight_layout()
    fig.savefig(output_dir / "pca_scatter_time.png", dpi=180)
    plt.close(fig)

    if "cluster" in score_table.columns and score_table["cluster"].notna().any():
        fig, axis = plt.subplots(figsize=(8, 6))
        for cluster_id, group in score_table.groupby("cluster", dropna=True):
            axis.scatter(group["PC1"], group["PC2"], s=8, alpha=0.45, label=str(cluster_id))
        axis.legend(title="Cluster", frameon=False, markerscale=2)
        axis.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
        axis.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
        axis.set_title("Per-cell PCA colored by k-means cluster")
        fig.tight_layout()
        fig.savefig(output_dir / "pca_scatter_cluster.png", dpi=180)
        plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 4))
    axis.plot(
        np.arange(1, n_components + 1),
        np.cumsum(pca.explained_variance_ratio_) * 100.0,
        marker="o",
    )
    axis.set_xlabel("Number of components")
    axis.set_ylabel("Cumulative explained variance (%)")
    axis.set_title("PCA scree (cumulative)")
    fig.tight_layout()
    fig.savefig(output_dir / "pca_scree.png", dpi=180)
    plt.close(fig)

    print(f"Cells: {len(table)}")
    print(f"Features used: {len(keep)}")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
