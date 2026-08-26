#!/usr/bin/env bash
#SBATCH --job-name=morph-rna-nulls
#SBATCH --partition=elscn.q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/morph-rna-nulls-%j.out
#SBATCH --error=logs/morph-rna-nulls-%j.err
#
# STUB: label-shuffle / null validation for morph→RNA program recovery.
# Not implemented yet — exits non-zero so accidental submits are obvious.
#
# Intended later:
#   - shuffle true_time (or P(time) rows) vs frozen atlas
#   - compare null mean Pearson to observed rna_metrics.json
#
#   N_NULLS=200 ROOT=... sbatch .../slurm/validation_nulls.slurm

set -euo pipefail

ROOT="${ROOT:-/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq}"
cd "$ROOT"
mkdir -p logs

echo "Host: $(hostname)  Job: ${SLURM_JOB_ID:-local}"
echo "ROOT=$ROOT  N_NULLS=${N_NULLS:-200}"
date

echo "ERROR: validation nulls not implemented yet." >&2
echo "Wire a scripts/run_validation_nulls.py (or flag on evaluate_rna_readout.py) before using this job." >&2
exit 1
