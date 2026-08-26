#!/bin/bash
set -euo pipefail
ROOT=/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq
RNA=/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2
mkdir -p "$ROOT"/{scripts,data,results,logs}
cd "$ROOT"

if [[ ! -x "$ROOT/.venv/bin/python" ]]; then
  python3 -m venv "$ROOT/.venv"
fi
"$ROOT/.venv/bin/pip" install --upgrade pip
"$ROOT/.venv/bin/pip" install "numpy<2.1" "scikit-learn" "matplotlib" "PyYAML" "joblib"

module load R/4.5.2 2>/dev/null || true
export R_LIBS_USER="${RNA}/R_library:${R_LIBS_USER:-}"
Rscript "$ROOT/scripts/export_vst.R" \
  "$RNA/data/vst_matrix.rds" \
  "$RNA/data/colData.rds" \
  "$ROOT/data/vst_matrix.csv" \
  "$ROOT/data/colData.csv"

"$ROOT/.venv/bin/python" "$ROOT/scripts/run_morph_to_rnaseq.py" --project-root "$ROOT"
