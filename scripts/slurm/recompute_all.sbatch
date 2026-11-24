#!/bin/bash
# Recompute every aggregate, statistic, and figure from cached per-spec runs,
# e.g. after a change to the S3 definition. Expects curated runs under
# results/cur_{hne,breast,fluo}/runs and factorial runs under results/full_*/runs.
#SBATCH --job-name=mosaic-recompute
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --partition=mit_normal

set -uo pipefail
source .venv/bin/activate
export PYTHONPATH=. PYTHONUNBUFFERED=1
step() { echo "=== $*"; "$@" || echo "!!! FAILED: $*"; }

# Point the default results folders (which every analysis script reads) at the
# curated runs.
link() { if [ -L "$2" ] || [ ! -e "$2" ]; then mkdir -p "$(dirname "$2")"; ln -sfn "$1" "$2"; else echo "!!! $2 exists and is not a link"; exit 1; fi; }
link "$PWD/results/cur_hne/runs" results/runs
link "$PWD/results/cur_breast/runs" results/breast_cancer/runs
link "$PWD/results/cur_fluo/runs" results/visium_fluo/runs

C=configs/multiverse_spec.yaml
step python scripts/main.py --config $C --lean
step python scripts/main.py --config $C --dataset breast_cancer --lean
step python scripts/main.py --config $C --dataset visium_fluo --lean
step python scripts/analyze_dataset.py --lean
step python scripts/analyze_dataset.py --dataset breast_cancer --lean
step python scripts/analyze_pathology.py
step python scripts/pathology_gaps.py
step python -m mosaic.consensus
step python -m mosaic.consensus_selection
step python -m mosaic.predict_fragility
step python -m mosaic.reference_robustness
step python -m mosaic.downstream
step python scripts/generate_figures.py
step python scripts/generate_tissue.py
step python scripts/generate_cross_dataset.py

F=configs/multiverse_full.yaml
step python scripts/main.py --config $F --results-root results/full_hne --lean
step python scripts/main.py --config $F --dataset breast_cancer --results-root results/full_breast_cancer --lean
step python scripts/analyze_dataset.py --results-root results/full_hne --lean
step python scripts/analyze_dataset.py --dataset breast_cancer --results-root results/full_breast_cancer --lean
step python scripts/variance_categorical.py results/aggregate/specification_curves.json \
    results/breast_cancer/aggregate/specification_curves.json \
    results/full_hne/aggregate/specification_curves.json \
    results/full_breast_cancer/aggregate/specification_curves.json

tar czf logs/recompute_outputs.tgz results/aggregate results/breast_cancer/aggregate \
    results/visium_fluo/aggregate results/full_hne/aggregate \
    results/full_breast_cancer/aggregate results/figures results/breast_cancer/figures \
    results/visium_fluo/figures 2>/dev/null
echo "=== done"
