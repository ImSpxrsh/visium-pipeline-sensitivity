"""
Determinism test for src/pipeline.py.

Marked ``network`` (downloads the Visium demo) and slow. Run with:
    pytest -m network
"""

import numpy as np
import pytest

SPEC = {
    "normalization": "log1p",
    "feature_selection": {"method": "seurat", "n_genes": 3000},
    "dimensionality_reduction": {"n_components": 30},
    "neighborhood_graph": {"type": "knn", "k_or_radius": 10},
    "clustering": {"method": "leiden", "resolution": 0.5},
    "spatially_variable_genes": {"method": "morans_i"},
}


@pytest.mark.network
def test_pipeline_determinism():
    from src.data_loader import load_dataset, prepare_data
    from src.pipeline import run_pipeline

    adata, _ = load_dataset("visium_hne", verbose=False)
    adata = prepare_data(adata, verbose=False)

    r1 = run_pipeline(adata.copy(), SPEC, random_state=42, verbose=False)
    r2 = run_pipeline(adata.copy(), SPEC, random_state=42, verbose=False)

    assert np.array_equal(
        r1["adata"].obs["leiden"].values, r2["adata"].obs["leiden"].values
    )
    assert r1["n_clusters"] == r2["n_clusters"]
