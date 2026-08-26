#!/usr/bin/env python3
"""Widescreen Design slide figure: morph → RNA pipeline (visual, few words)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parents[1] / "slides" / "assets" / "generated" / "design_pipeline.png"


def box(ax, x, y, w, h, text, fc="#f4f1ea", ec="#2c2c2c", ts=11, bold=False):
    p = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.4,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(p)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=ts,
        fontweight="bold" if bold else "normal",
        color="#1a1a1a",
        wrap=True,
    )


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.6,
            color="#333333",
        )
    )


def main() -> None:
    fig, ax = plt.subplots(figsize=(13.333, 7.5), dpi=160)
    ax.set_xlim(0, 13.333)
    ax.set_ylim(0, 7.5)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Title
    ax.text(0.45, 7.05, "Design: morphology  →  RNA", fontsize=26, fontweight="bold", color="#1a1a1a", va="center")
    ax.text(
        0.45,
        6.55,
        "Imaging predicts time · frozen bulk atlas supplies programs · no cell-matched RNA",
        fontsize=13,
        color="#555555",
        va="center",
    )

    # --- Row A: imaging prep ---
    ax.text(0.45, 5.95, "A. Imaging → slice features", fontsize=14, fontweight="bold", color="#333333")
    bw, bh = 2.7, 0.85
    yA = 4.85
    xs = [0.45, 3.55, 6.65, 9.75]
    labels_A = [
        "Cellpose\nsegmentation",
        "Single-cell\ncrops",
        "CellProfiler +\nskeleton metrics",
        "Pool to slice\n(median / p90 / clusters)",
    ]
    for x, lab in zip(xs, labels_A):
        box(ax, x, yA, bw, bh, lab, fc="#eef2f7")
    for i in range(3):
        arrow(ax, xs[i] + bw, yA + bh / 2, xs[i + 1], yA + bh / 2)

    # --- Row B: two-stage model ---
    ax.text(0.45, 4.35, "B. Two-stage model (inference uses morph only)", fontsize=14, fontweight="bold", color="#333333")
    yB = 3.15
    box(ax, 0.45, yB, 2.9, 0.95, "Slice morphology\nfeatures", fc="#e8f0e8", bold=True)
    arrow(ax, 3.35, yB + 0.48, 3.85, yB + 0.48)
    box(ax, 3.85, yB, 3.4, 0.95, "Stage 1\nMorph → P(time)\nL2 logistic · nested CV ~92%", fc="#f7ebe0", bold=True)
    arrow(ax, 7.25, yB + 0.48, 7.75, yB + 0.48)
    box(ax, 7.75, yB, 2.4, 0.95, "Mix\nP(time) · atlas", fc="#efe6f5", bold=True)
    arrow(ax, 10.15, yB + 0.48, 10.65, yB + 0.48)
    box(ax, 10.65, yB, 2.2, 0.95, "Predicted RNA\nprograms + genes", fc="#e8f0e8", bold=True)

    # --- Side / bottom: atlas built once ---
    ax.text(0.45, 2.65, "C. Frozen RNA atlas (built once from bulk RNA-seq)", fontsize=14, fontweight="bold", color="#333333")
    yC = 1.35
    box(ax, 0.45, yC, 3.2, 1.05, "Bulk wells\nDESeq2 → VST", fc="#f5f0e6")
    arrow(ax, 3.65, yC + 0.52, 4.15, yC + 0.52)
    box(ax, 4.15, yC, 4.0, 1.05, "Gene-set programs\n(+ panel genes)\nz-scored by time", fc="#f5f0e6")
    arrow(ax, 8.15, yC + 0.52, 8.65, yC + 0.52)
    box(ax, 8.65, yC, 4.2, 1.05, "Atlas lookup table\n5 times × K programs\nNot used as input on new FOVs", fc="#f5f0e6")

    # Dotted link from atlas up into mix
    ax.annotate(
        "",
        xy=(8.95, 3.15),
        xytext=(10.0, 2.4),
        arrowprops=dict(arrowstyle="-|>", color="#7a5a9a", lw=1.5, ls="--"),
    )
    ax.text(9.55, 2.85, "feeds Stage 2", fontsize=10, color="#7a5a9a", style="italic")

    # Footer constraints
    ax.text(
        0.45,
        0.55,
        "Constraints: imaging slices & bulk wells share time labels only  ·  no cell-matched transcriptomes  ·  new FOVs → morph CSV only",
        fontsize=11,
        color="#555555",
        va="center",
    )
    ax.text(
        0.45,
        0.2,
        "Replace the text-heavy Design slide with this figure (full bleed).",
        fontsize=9,
        color="#999999",
        va="center",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
