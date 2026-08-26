#!/usr/bin/env bash
#SBATCH --job-name=astro-segment
#SBATCH --partition=gpu.q
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# output/error overridden by submit_pipeline.sh → ROOT/logs/slurm/

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${ROOT:-$SLURM_SUBMIT_DIR}}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc_unified.yaml}"
export ROOT="${ROOT:-$REPO_ROOT}"

cd "$REPO_ROOT"
module load python/3.11.13
source "${REPO_ROOT}/envs/.venv/bin/activate" 2>/dev/null || source "${REPO_ROOT}/.venv/bin/activate"

astrocyte-pipeline inventory --config "$CONFIG"
astrocyte-pipeline segment --config "$CONFIG"
astrocyte-pipeline crops --config "$CONFIG"
