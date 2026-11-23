#!/bin/bash
# Statistics (variance decomposition, stakes, confound, convergence) and figures
# for a finished sweep. Run after aggregate.sbatch, e.g.:
#   sbatch --dependency=afterok:<AGG_JOB_ID> scripts/slurm/analyze.sbatch visium_hne results/full_hne
#SBATCH --job-name=mosaic-analyze
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --partition=mit_normal

set -euo pipefail
source .venv/bin/activate
export PYTHONUNBUFFERED=1
PYTHONPATH=. python scripts/analyze_dataset.py --dataset "$1" --results-root "$2" --lean
