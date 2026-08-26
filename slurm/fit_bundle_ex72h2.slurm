#!/usr/bin/env bash
#SBATCH --job-name=fit-bundle-ex72h2
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/logs/slurm/%x-%j.out
#SBATCH --error=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/logs/slurm/%x-%j.err

set -euo pipefail

ROOT=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end
TAL=/ems/elsc-labs/habib-n/tal.kupferman/morph_to_rnaseq_ex72h2
NESTED_SRC="$TAL/results_nested_top_features_ex72h2"
NESTED_DST="$ROOT/models/morph_to_rnaseq/nested"
BUNDLE="$ROOT/models/morph_to_rnaseq/bundle"
MORPH_CSV="$ROOT/models/morph_to_rnaseq/single_cell_full_features.csv"
# Prefer ROOT scripts; fall back to segev morph_to_rnaseq scripts
SCRIPTS="$ROOT/scripts/analysis/morph_to_rnaseq"
if [[ ! -f "$SCRIPTS/predict_rna_from_morph.py" ]]; then
  SCRIPTS=/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/scripts
fi
VENV=/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/.venv
if [[ ! -x "$VENV/bin/python" ]]; then
  VENV="$ROOT/envs/.venv"
fi

module load python/3.11.13 2>/dev/null || true
mkdir -p "$ROOT/logs/slurm" "$NESTED_DST" "$BUNDLE"

echo "=== copy nested results from tal ex72h2 ==="
rm -rf "$NESTED_DST"
mkdir -p "$NESTED_DST"
cp -a "$NESTED_SRC"/. "$NESTED_DST/"

echo "=== rebuild rna_by_time.npz from ex72h2 VST (if evaluate available) ==="
source "$VENV/bin/activate"
# Stage VST/colData names evaluate expects under a temp project layout
TMP_PROJ="$ROOT/models/morph_to_rnaseq/_fit_ex72h2_proj"
mkdir -p "$TMP_PROJ/data" "$TMP_PROJ/scripts" "$TMP_PROJ/results_nested_top_features"
cp -a "$TAL/data/vst_matrix_ex72h2.csv" "$TMP_PROJ/data/vst_matrix.csv"
cp -a "$TAL/data/colData_ex72h2.csv" "$TMP_PROJ/data/colData.csv"
cp -a "$NESTED_DST/programs.yaml" "$TMP_PROJ/scripts/programs.yaml"
# Point nested feature table for fit + evaluate (needs cv_predictions)
cp -a "$NESTED_DST/." "$TMP_PROJ/results_nested_top_features/"

# evaluate_rna_readout writes rna_by_time.npz into results_nested_top_features
if [[ -f "$SCRIPTS/evaluate_rna_readout.py" ]]; then
  python "$SCRIPTS/evaluate_rna_readout.py" \
    --project-root "$TMP_PROJ" \
    --nested-dir "$TMP_PROJ/results_nested_top_features" \
    --programs-yaml "$TMP_PROJ/scripts/programs.yaml" \
    --vst-csv "$TMP_PROJ/data/vst_matrix.csv" \
    --coldata-csv "$TMP_PROJ/data/colData.csv"
fi

# If evaluate wrote npz into TMP nested, copy to NESTED_DST and BUNDLE prep
if [[ -f "$TMP_PROJ/results_nested_top_features/rna_by_time.npz" ]]; then
  cp -a "$TMP_PROJ/results_nested_top_features/rna_by_time.npz" "$NESTED_DST/rna_by_time.npz"
  cp -a "$TMP_PROJ/results_nested_top_features/rna_metrics.json" "$NESTED_DST/" 2>/dev/null || true
elif [[ -f "$NESTED_DST/rna_by_time.npz" ]]; then
  echo "using existing nested rna_by_time.npz"
else
  # fallback: previous ROOT atlas
  echo "WARN: evaluate did not produce rna_by_time.npz; using prior ROOT bundle atlas if present"
  cp -a "$BUNDLE/rna_by_time.npz" "$NESTED_DST/rna_by_time.npz"
fi

# Fit expects nested dir layout with nested_feature_frequency + rna_by_time.npz
cp -a "$NESTED_DST/nested_feature_frequency.csv" "$TMP_PROJ/results_nested_top_features/"
cp -a "$NESTED_DST/rna_by_time.npz" "$TMP_PROJ/results_nested_top_features/"

echo "=== fit frozen bundle ==="
python "$SCRIPTS/predict_rna_from_morph.py" fit \
  --project-root "$TMP_PROJ" \
  --morph-csv "$MORPH_CSV" \
  --nested-dir "$TMP_PROJ/results_nested_top_features" \
  --bundle-dir "$BUNDLE"

echo "=== done ==="
ls -lh "$BUNDLE"
ls -lh "$NESTED_DST/nested_feature_frequency.csv" "$NESTED_DST/metrics.json" "$NESTED_DST/rna_by_time.npz"
python - <<'PY'
import json
from pathlib import Path
p = Path("/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/models/morph_to_rnaseq/nested/metrics.json")
print(json.dumps(json.loads(p.read_text()), indent=2)[:500])
PY
