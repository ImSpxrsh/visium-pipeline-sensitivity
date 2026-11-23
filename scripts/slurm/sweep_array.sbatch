#!/bin/bash
# Each array task runs PER_TASK consecutive specs (default 8), which keeps the
# job count under the cluster's per-user submit limit (448 on mit_normal).
# 864 specs / 8 = 108 tasks. Submit from the repo root, e.g.:
#   sbatch --array=0-107 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
# Args: <config> <dataset> <results-root> [per-task]
#SBATCH --job-name=mosaic-spec
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/%x_%A_%a.out
# Partition on MIT ORCD Engaging (12 h limit); change if `sinfo` lists others:
#SBATCH --partition=mit_normal

set -euo pipefail
CONFIG=${1:?config path}
DATASET=${2:?dataset name}
OUT=${3:?results root}
PER_TASK=${4:-8}

source .venv/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1   # write log lines as they happen
N_SPECS=$(grep -cE "^[[:space:]]*- name:" "$CONFIG")
START=$((SLURM_ARRAY_TASK_ID * PER_TASK))
END=$((START + PER_TASK - 1))
if [ "$END" -ge "$N_SPECS" ]; then END=$((N_SPECS - 1)); fi
for i in $(seq "$START" "$END"); do
    # Already-finished specs are skipped by content hash, so resubmitting is safe.
    python scripts/main.py --config "$CONFIG" --dataset "$DATASET" \
        --results-root "$OUT" --spec-index "$i" || echo "spec $i exited with an error"
done
