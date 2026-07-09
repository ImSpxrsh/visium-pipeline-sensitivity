"""
Verify the pipeline is deterministic: same seed + same spec -> identical clusters.

Run from the project root:  python scripts/test_reproducibility.py
Exits 0 if all runs agree, 1 otherwise.
"""

import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.data_loader import load_dataset, prepare_data  # noqa: E402
from mosaic.pipeline import run_pipeline  # noqa: E402

logger = logging.getLogger(__name__)

SPEC = {
    "normalization": "log1p",
    "feature_selection": {"method": "seurat", "n_genes": 3000},
    "dimensionality_reduction": {"n_components": 30},
    "neighborhood_graph": {"type": "knn", "k_or_radius": 10},
    "clustering": {"method": "leiden", "resolution": 0.5},
    "spatially_variable_genes": {"method": "morans_i"},
}


def test_reproducibility(n_runs: int = 3) -> bool:
    logging.basicConfig(level=logging.INFO)
    logger.info("Testing reproducibility with %d runs...", n_runs)

    adata, _ = load_dataset("visium_hne", verbose=False)
    adata = prepare_data(adata, verbose=False)

    cluster_runs = []
    for run_id in range(n_runs):
        logger.info("Run %d/%d ...", run_id + 1, n_runs)
        result = run_pipeline(adata.copy(), SPEC, random_state=42, verbose=False)
        cluster_runs.append(result["adata"].obs["leiden"].values)

    all_identical = True
    for i in range(1, len(cluster_runs)):
        if np.array_equal(cluster_runs[0], cluster_runs[i]):
            logger.info("Run 0 == Run %d: identical", i)
        else:
            logger.error("Run 0 != Run %d: clusters differ", i)
            all_identical = False

    logger.info("Reproducible" if all_identical else "NOT reproducible")
    return all_identical


if __name__ == "__main__":
    sys.exit(0 if test_reproducibility(n_runs=3) else 1)
