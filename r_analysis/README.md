# R morphology analysis

Downstream R analysis of CellProfiler / single-cell astrocyte feature tables. This folder is **separate from** the Python feature-extraction pipeline; it assumes feature CSVs already exist.

## Contents

| File | Purpose |
|------|---------|
| `astro_analysis_helpers.R` | Shared helpers (load, trimmed means, heatmaps, Kruskal–Wallis ranking) |
| `CP_AstrocyteMorphologyAnalysis.Rmd` | All five time points (`ctrl`, `4h`, `24h`, `72h`, `7d`) |
| `CP_AstrocyteMorphology_3Conditions.Rmd` | Focused `ctrl` / `24h` / `72h` comparison |
| `CP_SingleCellFeaturesHeatmap.R` | Heatmap for `single_cell_features_and_clusters.csv` |
| `CP_IntensityBoxplots.R` | Intensity-only Time boxplots from `AstroResultsAstrocytes_with_nuclei.csv` |

## Requirements

R packages: `tidyverse`, `pheatmap`. Pandoc is optional (notebooks can be sourced chunk-by-chunk with `Rscript`).

## Paths

Edit `DATA_PATH` and `OUTPUT_ROOT` (or `OUTPUT_DIR`) at the top of each script / Rmd to match your machine before running.

Typical inputs:

- CellProfiler object table (e.g. `AstroResultsAstrocytes.csv`)
- Optional: `single_cell_features_and_clusters.csv` for the standalone heatmap script

## Outputs

Plots are written under the configured output directory (heatmaps, stats, boxplots). No CSVs are required for the kept figures.
