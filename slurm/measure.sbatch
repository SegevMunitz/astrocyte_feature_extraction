#!/usr/bin/env bash
#SBATCH --job-name=astro-measure
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$SLURM_SUBMIT_DIR}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/elsc.yaml}"

cd "$REPO_ROOT"
module load python/3.11.13 java/21.0.9
source .venv/bin/activate

astrocyte-pipeline measure --config "$CONFIG"
