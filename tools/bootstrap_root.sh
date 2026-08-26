#!/usr/bin/env bash
# Curated bootstrap for ROOT=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end
# Copies only allowlisted runtime assets. Does NOT copy images.
# Does NOT submit jobs or run Drive FOVs.
set -euo pipefail

ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end}"
REPO="${REPO:-/ems/elsc-labs/habib-n/segev.munitz/astrocyte_feature_extraction}"
SRC_CKPT="/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/outputs/cellpose/three_channel/three_channel_transfer_20260816/lr_0p003/result/models/cellpose_3ch_lr_0p003_epoch_0250"
SRC_BUNDLE="/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq/bundle"
SRC_CLUSTER_CSV="/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results/single_cell_features_and_clusters.csv"
SRC_MORPH_CSV="/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results/single_cell_full_features.csv"
MANIFEST="$ROOT/logs/bootstrap_manifest.txt"

mkdir -p \
  "$ROOT/models/cellpose" \
  "$ROOT/models/morph_to_rnaseq/bundle" \
  "$ROOT/reference" \
  "$ROOT/cellprofiler" \
  "$ROOT/configs" \
  "$ROOT/scripts/analysis/morph_to_rnaseq" \
  "$ROOT/slurm" \
  "$ROOT/envs" \
  "$ROOT/outputs/runs" \
  "$ROOT/logs/slurm" \
  "$ROOT/tools"

copy_one() {
  local src="$1" dst="$2"
  if [[ ! -e "$src" ]]; then
    echo "MISSING: $src" | tee -a "$MANIFEST"
    return 1
  fi
  mkdir -p "$(dirname "$dst")"
  cp -a "$src" "$dst"
  local bytes
  bytes=$(stat -c%s "$dst" 2>/dev/null || stat -f%z "$dst")
  echo "COPIED $bytes $src -> $dst" | tee -a "$MANIFEST"
}

echo "=== bootstrap ROOT=$ROOT ===" | tee "$MANIFEST"
# --- Cellpose: single checkpoint only ---
copy_one "$SRC_CKPT" "$ROOT/models/cellpose/cellpose_3ch_lr_0p003_epoch_0250"

# --- Morph→RNA bundle allowlist ---
copy_one "$SRC_BUNDLE/model_bundle.joblib" "$ROOT/models/morph_to_rnaseq/bundle/model_bundle.joblib"
copy_one "$SRC_BUNDLE/rna_by_time.npz" "$ROOT/models/morph_to_rnaseq/bundle/rna_by_time.npz"
if [[ -f "$SRC_BUNDLE/README.txt" ]]; then
  copy_one "$SRC_BUNDLE/README.txt" "$ROOT/models/morph_to_rnaseq/bundle/README.txt"
fi

# --- Cluster training table (for assigner + optional join) ---
if [[ -f "$SRC_CLUSTER_CSV" ]]; then
  copy_one "$SRC_CLUSTER_CSV" "$ROOT/models/morph_to_rnaseq/single_cell_features_and_clusters.csv"
else
  echo "WARN: cluster CSV missing at $SRC_CLUSTER_CSV" | tee -a "$MANIFEST"
fi

# --- Reference / CP pipeline ---
copy_one "$REPO/reference/AstroResultsAstrocytes_with_nuclei.csv" \
  "$ROOT/reference/AstroResultsAstrocytes_with_nuclei.csv" || true
# reference dir may only exist on cluster repo
if [[ -d "$REPO/reference" ]]; then
  # copy only csv + cppipe if present
  find "$REPO/reference" -maxdepth 1 \( -name '*.csv' -o -name '*.cppipe' \) -type f | while read -r f; do
    copy_one "$f" "$ROOT/reference/$(basename "$f")"
  done
fi
if [[ -f "$REPO/cellprofiler/Astrocytes_CellposeLabels.cppipe" ]]; then
  copy_one "$REPO/cellprofiler/Astrocytes_CellposeLabels.cppipe" \
    "$ROOT/cellprofiler/Astrocytes_CellposeLabels.cppipe"
fi

# --- Code sync (allowlisted trees) ---
rsync -a --delete \
  --exclude '__pycache__' --exclude '*.pyc' --exclude '.pytest_cache' \
  "$REPO/src/" "$ROOT/src/"
rsync -a \
  --exclude '__pycache__' --exclude '*.pyc' \
  "$REPO/scripts/analysis/measure_single_cell_crops_full.py" "$ROOT/scripts/analysis/"
rsync -a \
  --exclude '__pycache__' --exclude '*.pyc' \
  "$REPO/scripts/analysis/morph_to_rnaseq/" "$ROOT/scripts/analysis/morph_to_rnaseq/"
rsync -a "$REPO/scripts/run_cellprofiler_direct.py" "$ROOT/scripts/" 2>/dev/null || true
rsync -a "$REPO/configs/elsc_unified.yaml" "$ROOT/configs/"
rsync -a "$REPO/slurm/" "$ROOT/slurm/"
rsync -a "$REPO/tools/bootstrap_root.sh" "$ROOT/tools/" 2>/dev/null || true
rsync -a "$REPO/pyproject.toml" "$ROOT/" 2>/dev/null || true

# --- Envs: prefer symlink/copy existing rather than redownload Cellpose from Tal ---
# Pipeline venv
if [[ ! -d "$ROOT/envs/.venv" ]]; then
  if [[ -d "$REPO/.venv" ]]; then
    echo "Linking pipeline venv from repo .venv" | tee -a "$MANIFEST"
    ln -sfn "$REPO/.venv" "$ROOT/envs/.venv"
  fi
fi
# CellProfiler venv
if [[ ! -d "$ROOT/envs/.cpvenv" ]]; then
  if [[ -d "$REPO/.cpvenv" ]]; then
    echo "Linking CP venv from repo .cpvenv" | tee -a "$MANIFEST"
    ln -sfn "$REPO/.cpvenv" "$ROOT/envs/.cpvenv"
  fi
fi

# Convenience top-level venv paths expected by some sbatch scripts
ln -sfn "$ROOT/envs/.venv" "$ROOT/.venv" 2>/dev/null || true
ln -sfn "$ROOT/envs/.cpvenv" "$ROOT/.cpvenv" 2>/dev/null || true

# --- Fit cluster assigner if morph+cluster labels available ---
module load python/3.11.13 2>/dev/null || true
if [[ -x "$ROOT/envs/.venv/bin/python" && -f "$SRC_MORPH_CSV" ]]; then
  "$ROOT/envs/.venv/bin/python" "$ROOT/scripts/analysis/morph_to_rnaseq/fit_cluster_assigner.py" \
    --morph-csv "$SRC_MORPH_CSV" \
    --out "$ROOT/models/morph_to_rnaseq/cluster_assigner.joblib" \
    | tee -a "$MANIFEST" || echo "WARN: cluster assigner fit failed" | tee -a "$MANIFEST"
elif [[ -x "$ROOT/envs/.venv/bin/python" && -f "$ROOT/models/morph_to_rnaseq/single_cell_features_and_clusters.csv" ]]; then
  "$ROOT/envs/.venv/bin/python" "$ROOT/scripts/analysis/morph_to_rnaseq/fit_cluster_assigner.py" \
    --morph-csv "$ROOT/models/morph_to_rnaseq/single_cell_features_and_clusters.csv" \
    --out "$ROOT/models/morph_to_rnaseq/cluster_assigner.joblib" \
    | tee -a "$MANIFEST" || echo "WARN: cluster assigner fit failed" | tee -a "$MANIFEST"
fi

# Editable install into ROOT venv if present
if [[ -x "$ROOT/envs/.venv/bin/pip" && -f "$ROOT/pyproject.toml" ]]; then
  "$ROOT/envs/.venv/bin/pip" install -e "$ROOT" -q || true
fi

echo "NOTE: images are NOT copied; config images_dir points at existing test_images." | tee -a "$MANIFEST"
echo "NOTE: Tal Cellpose env may still be used via linked/shared venv if torch/cellpose live there." | tee -a "$MANIFEST"
echo "DONE. Manifest: $MANIFEST" | tee -a "$MANIFEST"
