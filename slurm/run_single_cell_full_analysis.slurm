#!/usr/bin/env bash
#SBATCH --job-name=single-cell-full
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/single-cell-full-%j.out

set -euo pipefail

ANALYSIS_ROOT="${ANALYSIS_ROOT:-/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis}"
REPO_ROOT="${REPO_ROOT:-/ems/elsc-labs/habib-n/segev.munitz/astrocyte_feature_extraction}"

cd "$ANALYSIS_ROOT"
mkdir -p logs

module load python/3.11.13 java/21.0.9

if [[ -f "$ANALYSIS_ROOT/.venv_astro_pipeline/bin/activate" ]]; then
  source "$ANALYSIS_ROOT/.venv_astro_pipeline/bin/activate"
elif [[ -f "$REPO_ROOT/.venv/bin/activate" ]]; then
  source "$REPO_ROOT/.venv/bin/activate"
else
  echo "No Python venv found under analysis or repo" >&2
  exit 1
fi

python "$REPO_ROOT/scripts/measure_single_cell_crops_full.py" \
  --analysis-root "$ANALYSIS_ROOT" \
  --existing-features "$ANALYSIS_ROOT/results/single_cell_features_and_clusters.csv"

python "$REPO_ROOT/scripts/single_cell_pca_analysis.py" \
  --input "$ANALYSIS_ROOT/results/single_cell_full_features.csv" \
  --output-dir "$ANALYSIS_ROOT/results/plots/per_cell_pca"

echo "Done. Full features: $ANALYSIS_ROOT/results/single_cell_full_features.csv"
echo "PCA plots: $ANALYSIS_ROOT/results/plots/per_cell_pca/"
