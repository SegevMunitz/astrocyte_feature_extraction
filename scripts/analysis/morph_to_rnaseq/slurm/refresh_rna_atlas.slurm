#!/usr/bin/env bash
#SBATCH --job-name=morph-rna-atlas
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/morph-rna-atlas-%j.out
#SBATCH --error=logs/morph-rna-atlas-%j.err
#
# Rebuild frozen RNA atlas (evaluate_rna_readout) and refit the morph→RNA bundle.
# Submit from the morph_to_rnaseq project root (or any cwd with a logs/ dir):
#   sbatch /path/to/repo/scripts/analysis/morph_to_rnaseq/slurm/refresh_rna_atlas.slurm
# Or after syncing scripts onto the cluster tree:
#   sbatch $ROOT/scripts/slurm/refresh_rna_atlas.slurm
#
# Optional env overrides:
#   ROOT=/ems/.../morph_to_rnaseq
#   RNA=/ems/.../051223_BMP4_deseq2
#   SKIP_FIT=1   # evaluate only

set -euo pipefail

ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq}"
RNA="${RNA:-/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2}"
PY="${ROOT}/.venv/bin/python"

cd "$ROOT"
mkdir -p logs

echo "Host: $(hostname)  Job: ${SLURM_JOB_ID:-local}"
echo "ROOT=$ROOT"
echo "RNA=$RNA"
date

if [[ ! -x "$PY" ]]; then
  echo "ERROR: missing venv python at $PY" >&2
  echo "Create once on a compute node or via setup, then re-submit." >&2
  exit 127
fi

"$PY" "$ROOT/scripts/evaluate_rna_readout.py" \
  --project-root "$ROOT" \
  --programs-yaml "$ROOT/scripts/programs.yaml" \
  --contrast-dir "$RNA/results/contrasts"

if [[ "${SKIP_FIT:-0}" != "1" ]]; then
  "$PY" "$ROOT/scripts/predict_rna_from_morph.py" fit \
    --project-root "$ROOT"
fi

echo "Refreshed RNA atlas (and bundle unless SKIP_FIT=1) under $ROOT"
date
