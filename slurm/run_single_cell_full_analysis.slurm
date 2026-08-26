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
SCRIPTS_DIR="${SCRIPTS_DIR:-$ANALYSIS_ROOT/scripts}"

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

if [[ ! -f "$SCRIPTS_DIR/measure_single_cell_crops_full.py" ]]; then
  if [[ -f "$REPO_ROOT/scripts/measure_single_cell_crops_full.py" ]]; then
    SCRIPTS_DIR="$REPO_ROOT/scripts"
  else
    echo "Could not find measure_single_cell_crops_full.py under:" >&2
    echo "  $ANALYSIS_ROOT/scripts" >&2
    echo "  $REPO_ROOT/scripts" >&2
    exit 1
  fi
fi

echo "Python: $(which python)"
python -V
python - <<'PY'
import importlib.util
import sys

for package in ("cellprofiler", "cellprofiler_core", "sklearn", "matplotlib"):
    spec = importlib.util.find_spec(package)
    print(f"{package}: {'OK' if spec else 'MISSING'}")
    if package.startswith("cellprofiler") and spec is None:
        sys.exit("CellProfiler is required for measurement")
PY

echo "Scripts dir: $SCRIPTS_DIR"

python "$SCRIPTS_DIR/measure_single_cell_crops_full.py" \
  --analysis-root "$ANALYSIS_ROOT" \
  --existing-features "$ANALYSIS_ROOT/results/single_cell_features_and_clusters.csv" \
  --dry-run

python "$SCRIPTS_DIR/measure_single_cell_crops_full.py" \
  --analysis-root "$ANALYSIS_ROOT" \
  --existing-features "$ANALYSIS_ROOT/results/single_cell_features_and_clusters.csv"

python "$SCRIPTS_DIR/single_cell_pca_analysis.py" \
  --input "$ANALYSIS_ROOT/results/single_cell_full_features.csv" \
  --output-dir "$ANALYSIS_ROOT/results/plots/per_cell_pca"

echo "Done. Full features: $ANALYSIS_ROOT/results/single_cell_full_features.csv"
echo "PCA plots: $ANALYSIS_ROOT/results/plots/per_cell_pca/"
