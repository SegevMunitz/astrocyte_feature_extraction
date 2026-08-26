#!/usr/bin/env bash
# Refresh ROOT with latest morph CSV, nested morph→RNA artifacts, and finish install.
set -euo pipefail

ROOT=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end
REPO=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_feature_extraction
MORPH_PROJ=/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq
ANALYSIS=/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results
LATEST_MORPH="$ANALYSIS/single_cell_full_features.csv"
LATEST_CLUSTER="$ANALYSIS/single_cell_features_and_clusters.csv"
NESTED="$MORPH_PROJ/results_nested_top_features"
BUNDLE="$MORPH_PROJ/bundle"
LOG="$ROOT/logs/refresh_latest_assets.txt"

mkdir -p "$ROOT/models/morph_to_rnaseq/bundle" \
  "$ROOT/models/morph_to_rnaseq/nested" \
  "$ROOT/reference" \
  "$ROOT/logs" \
  "$ROOT/docs"

{
  echo "=== refresh $(date -Iseconds) ==="
  echo "morph: $LATEST_MORPH ($(stat -c%s "$LATEST_MORPH") bytes, $(stat -c%y "$LATEST_MORPH"))"
  echo "cluster: $LATEST_CLUSTER"
  echo "bundle: $BUNDLE"
  echo "nested: $NESTED"
} | tee "$LOG"

# Latest full cell feature matrix (Aug 23 run with CP+skeleton)
cp -a "$LATEST_MORPH" "$ROOT/models/morph_to_rnaseq/single_cell_full_features.csv"
cp -a "$LATEST_CLUSTER" "$ROOT/models/morph_to_rnaseq/single_cell_features_and_clusters.csv"

# Latest frozen bundle + nested CV sidecars used to build it
cp -a "$BUNDLE/model_bundle.joblib" "$BUNDLE/rna_by_time.npz" "$BUNDLE/README.txt" \
  "$ROOT/models/morph_to_rnaseq/bundle/"
for f in nested_feature_frequency.csv nested_fold_features.csv programs.yaml \
         rna_by_time.npz rna_metrics.json program_viability.json metrics.json \
         image_features_used.txt cell_features_used.txt; do
  if [[ -f "$NESTED/$f" ]]; then
    cp -a "$NESTED/$f" "$ROOT/models/morph_to_rnaseq/nested/$f"
    echo "COPIED nested/$f" | tee -a "$LOG"
  fi
done

# Reference schema: derive a lightweight header CSV from the latest morph matrix
# (FOV CP schema oracle can stay; also store morph header for crop-measure parity)
module load python/3.11.13
source "$ROOT/envs/.venv/bin/activate"
python - <<'PY'
from pathlib import Path
import csv
src = Path("/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/models/morph_to_rnaseq/single_cell_full_features.csv")
dst = Path("/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/reference/single_cell_full_features_header.csv")
with src.open(newline="", encoding="utf-8") as handle:
    cols = next(csv.reader(handle))
# Write a 1-row schema CSV (header only + dummy) for documentation / optional use
dst.parent.mkdir(parents=True, exist_ok=True)
with dst.open("w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle)
    w.writerow(cols)
    w.writerow([""] * len(cols))
print({"n_cols": len(cols), "has_cluster": "cluster" in cols, "out": str(dst)})
PY

# Keep AstroResults for FOV measure path, but also sync any newer cppipe from REPO
if [[ -f "$REPO/cellprofiler/Astrocytes_CellposeLabels.cppipe" ]]; then
  cp -a "$REPO/cellprofiler/Astrocytes_CellposeLabels.cppipe" "$ROOT/cellprofiler/"
fi
# Prefer morph project's programs.yaml into scripts if present
if [[ -f "$NESTED/programs.yaml" ]]; then
  mkdir -p "$ROOT/scripts/analysis/morph_to_rnaseq"
  cp -a "$NESTED/programs.yaml" "$ROOT/scripts/analysis/morph_to_rnaseq/programs.yaml"
fi

# Finish deps + cluster assigner on LATEST full features
pip install -q joblib scikit-learn
python "$ROOT/scripts/analysis/morph_to_rnaseq/fit_cluster_assigner.py" \
  --morph-csv "$ROOT/models/morph_to_rnaseq/single_cell_full_features.csv" \
  --out "$ROOT/models/morph_to_rnaseq/cluster_assigner.joblib" | tee -a "$LOG"

# Editable install
cp -a "$REPO/README.md" "$ROOT/README.md" 2>/dev/null || true
cp -a "$REPO/docs/UNIFIED_PIPELINE.md" "$ROOT/docs/" 2>/dev/null || true
pip install -e "$ROOT" -q || pip install -e "$REPO" -q

which astrocyte-pipeline || true
astrocyte-pipeline --help | head -25 || true
ls -la "$ROOT/models/morph_to_rnaseq/"
ls -la "$ROOT/models/morph_to_rnaseq/bundle/"
ls -la "$ROOT/models/morph_to_rnaseq/cluster_assigner.joblib"
ls -la "$ROOT/reference/"
echo REFRESH_OK | tee -a "$LOG"
