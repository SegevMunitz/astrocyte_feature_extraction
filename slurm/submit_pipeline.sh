#!/usr/bin/env bash
# Submit unified chain: segment(+crops) → crop_measure → predict_rna
# All Slurm logs under $ROOT/logs/slurm/
set -euo pipefail

ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end}"
CONFIG="${CONFIG:-$ROOT/configs/elsc_unified.yaml}"
LOGDIR="${LOGDIR:-$ROOT/logs/slurm}"
mkdir -p "$LOGDIR"

segment_job="$(
  sbatch --parsable \
    --output="$LOGDIR/%x-%j.out" \
    --error="$LOGDIR/%x-%j.err" \
    --export="ALL,REPO_ROOT=$ROOT,CONFIG=$CONFIG,ROOT=$ROOT" \
    "$ROOT/slurm/segment.sbatch"
)"
crop_job="$(
  sbatch --parsable \
    --dependency="afterok:$segment_job" \
    --output="$LOGDIR/%x-%j.out" \
    --error="$LOGDIR/%x-%j.err" \
    --export="ALL,REPO_ROOT=$ROOT,CONFIG=$CONFIG,ROOT=$ROOT" \
    "$ROOT/slurm/crop_measure.sbatch"
)"
predict_job="$(
  sbatch --parsable \
    --dependency="afterok:$crop_job" \
    --output="$LOGDIR/%x-%j.out" \
    --error="$LOGDIR/%x-%j.err" \
    --export="ALL,REPO_ROOT=$ROOT,CONFIG=$CONFIG,ROOT=$ROOT" \
    "$ROOT/slurm/predict_rna.sbatch"
)"

printf 'ROOT: %s\nCONFIG: %s\nLOGDIR: %s\n' "$ROOT" "$CONFIG" "$LOGDIR"
printf 'Segmentation(+crops) job: %s\nCrop-measure job: %s\nPredict-RNA job: %s\n' \
  "$segment_job" "$crop_job" "$predict_job"
