#!/bin/bash
# Rebuild the frozen RNA atlas from programs.yaml + existing nested CV.
# Does not retrain morph → time.
set -euo pipefail
ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq}"
RNA="${RNA:-/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2}"
HERE="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$ROOT/scripts" "$ROOT/data" "$ROOT/logs"
if [[ "$(cd "$HERE" && pwd)" != "$(cd "$ROOT/scripts" && pwd)" ]]; then
  cp -f "$HERE"/*.py "$HERE"/*.yaml "$HERE"/*.md "$ROOT/scripts/"
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
