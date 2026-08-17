#!/usr/bin/env bash
#SBATCH --job-name=astro-segment
#SBATCH --partition=gpu.q
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$SLURM_SUBMIT_DIR}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc.yaml}"

cd "$REPO_ROOT"
module load python/3.11.13
source .venv/bin/activate

astrocyte-pipeline inventory --config "$CONFIG"
astrocyte-pipeline segment --config "$CONFIG"
astrocyte-pipeline crops --config "$CONFIG"
astrocyte-pipeline spatial --config "$CONFIG"
