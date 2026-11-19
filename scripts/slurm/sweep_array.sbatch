#!/bin/bash
# One spec per job. Submit from the repo root, e.g.:
#   sbatch --array=0-863%50 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
# Args: <config> <dataset> <results-root>
#SBATCH --job-name=mosaic-spec
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm/%x_%A_%a.out
# Uncomment and set to a partition from `sinfo` if your cluster needs one:
##SBATCH --partition=mit_normal

set -euo pipefail
CONFIG=${1:?config path}
DATASET=${2:?dataset name}
OUT=${3:?results root}

source .venv/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python scripts/main.py --config "$CONFIG" --dataset "$DATASET" \
    --results-root "$OUT" --spec-index "$SLURM_ARRAY_TASK_ID"
