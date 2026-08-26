#!/usr/bin/env bash
#SBATCH --job-name=astro-predict-rna
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${ROOT:-$SLURM_SUBMIT_DIR}}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc_unified.yaml}"
export ROOT="${ROOT:-$REPO_ROOT}"

cd "$REPO_ROOT"
module load python/3.11.13
source "${REPO_ROOT}/envs/.venv/bin/activate" 2>/dev/null || source "${REPO_ROOT}/.venv/bin/activate"

astrocyte-pipeline predict_rna --config "$CONFIG"
