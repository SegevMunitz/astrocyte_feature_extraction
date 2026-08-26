#!/bin/bash
# Rebuild the frozen RNA atlas from programs.yaml + existing nested CV.
# Does not retrain morph → time.
#
# On the ELSC cluster: do NOT run this interactively on the login node.
# Submit via Slurm instead:
#   cd /ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq
#   sbatch scripts/slurm/refresh_rna_atlas.slurm
# (Repo path: scripts/analysis/morph_to_rnaseq/slurm/refresh_rna_atlas.slurm —
#  sync that file into $ROOT/scripts/slurm/ if needed.)
set -euo pipefail
ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq}"
RNA="${RNA:-/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2}"
HERE="$(cd "$(dirname "$0")" && pwd)"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "WARNING: not inside a Slurm allocation. On ELSC, prefer:" >&2
  echo "  sbatch scripts/slurm/refresh_rna_atlas.slurm" >&2
  echo "Continuing anyway (local/dev use)." >&2
fi

mkdir -p "$ROOT/scripts" "$ROOT/data" "$ROOT/logs" "$ROOT/scripts/slurm"
if [[ "$(cd "$HERE" && pwd)" != "$(cd "$ROOT/scripts" && pwd)" ]]; then
  cp -f "$HERE"/*.py "$HERE"/*.yaml "$HERE"/*.md "$ROOT/scripts/"
  if [[ -d "$HERE/slurm" ]]; then
    cp -f "$HERE"/slurm/*.slurm "$ROOT/scripts/slurm/"
  fi
fi

if [[ ! -x "$ROOT/.venv/bin/python" ]]; then
  python3 -m venv "$ROOT/.venv"
fi
"$ROOT/.venv/bin/pip" install -q "numpy<2.1" scikit-learn matplotlib PyYAML joblib

"$ROOT/.venv/bin/python" "$ROOT/scripts/evaluate_rna_readout.py" \
  --project-root "$ROOT" \
  --programs-yaml "$ROOT/scripts/programs.yaml" \
  --contrast-dir "$RNA/results/contrasts"

"$ROOT/.venv/bin/python" "$ROOT/scripts/predict_rna_from_morph.py" fit \
  --project-root "$ROOT"

echo "Refreshed RNA atlas and bundle under $ROOT"
