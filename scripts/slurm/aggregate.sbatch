#!/bin/bash
# Run after the array finishes. Loads every cached spec and computes S3, ARI, etc.
#   sbatch --dependency=afterany:<ARRAY_JOB_ID> scripts/slurm/aggregate.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
#SBATCH --job-name=mosaic-agg
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --partition=mit_normal

set -euo pipefail
source .venv/bin/activate
python scripts/main.py --config "$1" --dataset "$2" --results-root "$3" --lean
