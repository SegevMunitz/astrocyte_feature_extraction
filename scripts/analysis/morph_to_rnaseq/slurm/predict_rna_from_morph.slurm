#!/usr/bin/env bash
#SBATCH --job-name=morph-rna-pred
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/morph-rna-pred-%j.out
#SBATCH --error=logs/morph-rna-pred-%j.err
#
# Morph → RNA predict (default) or fit bundle.
# Submit from a directory that can hold logs/ (usually $ROOT):
#
#   # Predict on a new morphology CSV
#   MORPH_CSV=/path/to/new_cells.csv \
#   OUT_DIR=/path/to/pred_out \
#   sbatch .../slurm/predict_rna_from_morph.slurm
#
#   # Fit only (rebuild bundle from nested CV + morph CSV)
#   MODE=fit sbatch .../slurm/predict_rna_from_morph.slurm
#
# Env:
#   MODE=predict|fit   (default: predict)
#   ROOT, BUNDLE_DIR, MORPH_CSV, OUT_DIR, MIN_FOLDS, C

set -euo pipefail

ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq}"
PY="${ROOT}/.venv/bin/python"
MODE="${MODE:-predict}"
BUNDLE_DIR="${BUNDLE_DIR:-$ROOT/bundle}"
MORPH_CSV="${MORPH_CSV:-/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results/single_cell_full_features.csv}"

cd "$ROOT"
mkdir -p logs

echo "Host: $(hostname)  Job: ${SLURM_JOB_ID:-local}"
echo "MODE=$MODE  ROOT=$ROOT"
echo "MORPH_CSV=$MORPH_CSV"
date

if [[ ! -x "$PY" ]]; then
  echo "ERROR: missing venv python at $PY" >&2
  exit 127
fi

case "$MODE" in
  predict)
    if [[ -z "${MORPH_CSV:-}" || ! -f "$MORPH_CSV" ]]; then
      echo "ERROR: set MORPH_CSV to an existing morphology CSV for predict" >&2
      exit 2
    fi
    PRED_ARGS=(predict --morph-csv "$MORPH_CSV" --bundle-dir "$BUNDLE_DIR")
    if [[ -n "${OUT_DIR:-}" ]]; then
      PRED_ARGS+=(--out-dir "$OUT_DIR")
    fi
    "$PY" "$ROOT/scripts/predict_rna_from_morph.py" "${PRED_ARGS[@]}"
    ;;
  fit)
    FIT_ARGS=(fit --project-root "$ROOT" --morph-csv "$MORPH_CSV" --bundle-dir "$BUNDLE_DIR")
    if [[ -n "${MIN_FOLDS:-}" ]]; then
      FIT_ARGS+=(--min-folds "$MIN_FOLDS")
    fi
    if [[ -n "${C:-}" ]]; then
      FIT_ARGS+=(--C "$C")
    fi
    "$PY" "$ROOT/scripts/predict_rna_from_morph.py" "${FIT_ARGS[@]}"
    ;;
  *)
    echo "ERROR: MODE must be predict or fit (got: $MODE)" >&2
    exit 2
    ;;
esac

echo "Done ($MODE)."
date
