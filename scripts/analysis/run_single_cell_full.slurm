#!/usr/bin/env bash
#SBATCH --job-name=single-cell-full
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/single-cell-full-%j.out
#SBATCH --error=logs/single-cell-full-%j.err

set -euo pipefail

ANALYSIS_ROOT="/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis"
SIF="/ems/elsc-labs/habib-n/tal.kupferman/Segmentation App/containers/habib_cp_cellpose_4.2.8_cpu.sif"
SCRIPT="$ANALYSIS_ROOT/scripts/measure_single_cell_crops_full.py"

cd "$ANALYSIS_ROOT"
mkdir -p logs

module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_BIN="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_BIN="$(command -v singularity)"
elif [[ -x /ems/elsc-labs/apps/singularity/4.3.2/bin/singularity ]]; then
  CONTAINER_BIN="/ems/elsc-labs/apps/singularity/4.3.2/bin/singularity"
else
  echo "ERROR: apptainer/singularity not found on $(hostname)" >&2
  exit 127
fi

EXTRA_ARGS=(--workers "${WORKERS:-8}" --skip-existing)
if [[ -n "${LIMIT:-}" ]]; then
  EXTRA_ARGS+=(--limit "$LIMIT")
fi

echo "Host: $(hostname)  Job: ${SLURM_JOB_ID:-local}"
echo "Container: $CONTAINER_BIN"
echo "SIF: $SIF"
echo "Script: $SCRIPT"
echo "Extra args: ${EXTRA_ARGS[*]}"
date

"$CONTAINER_BIN" exec \
  --bind "$ANALYSIS_ROOT:$ANALYSIS_ROOT" \
  --bind "/ems/elsc-labs/habib-n/tal.kupferman:/ems/elsc-labs/habib-n/tal.kupferman" \
  "$SIF" \
  python "$SCRIPT" \
  --analysis-root "$ANALYSIS_ROOT" \
  --existing-features "$ANALYSIS_ROOT/results/single_cell_features_and_clusters.csv" \
  "${EXTRA_ARGS[@]}"

echo "Done."
date
