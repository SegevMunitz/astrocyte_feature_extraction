#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc.yaml}"

segment_job="$(
  sbatch --parsable \
    --export="ALL,REPO_ROOT=$REPO_ROOT,CONFIG=$CONFIG" \
    "$REPO_ROOT/slurm/segment.sbatch"
)"
measurement_job="$(
  sbatch --parsable \
    --dependency="afterok:$segment_job" \
    --export="ALL,REPO_ROOT=$REPO_ROOT,CONFIG=$CONFIG" \
    "$REPO_ROOT/slurm/measure.sbatch"
)"

printf 'Segmentation job: %s\nMeasurement job: %s\n' "$segment_job" "$measurement_job"
