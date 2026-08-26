#!/usr/bin/env bash
#SBATCH --job-name=astro-crop-measure
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${ROOT:-$SLURM_SUBMIT_DIR}}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc_unified.yaml}"
export ROOT="${ROOT:-$REPO_ROOT}"

cd "$REPO_ROOT"
module load python/3.11.13 java/21.0.9
source "${REPO_ROOT}/envs/.venv/bin/activate" 2>/dev/null || source "${REPO_ROOT}/.venv/bin/activate"

astrocyte-pipeline crop_measure --config "$CONFIG"
