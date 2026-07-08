"""
Unit tests for src/statistics.py on small synthetic sweeps (no network).
"""

import anndata as ad
import numpy as np
import pandas as pd

from src.statistics import (
    _s3_from_matrix,
    build_aligned_matrix,
    confound_check_umi,
    multiverse_size_convergence,
    variance_decomposition,
)


def _make_run(labels, spec, n_domains, regions=None, counts=None, n_genes=3):
    n = len(labels)
    a = ad.AnnData(np.ones((n, n_genes)))
    a.obs_names = [f"spot{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["leiden"] = pd.Categorical([str(x) for x in labels])
    if regions is not None:
        a.obs["cluster"] = pd.Categorical(regions)
    if counts is not None:
        a.obs["total_counts"] = counts
    a.obsm["spatial"] = np.random.RandomState(0).rand(n, 2)
    return {"adata": a, "spec": spec, "metadata": {"n_clusters_found": n_domains}}


def test_s3_from_matrix():
    mat = np.array([[0, 0, 0], [0, 1, 1]])  # spot0 unanimous, spot1 split 1/3 vs 2/3
    s3 = _s3_from_matrix(mat)
    assert s3[0] == 1.0
    assert abs(s3[1] - 2 / 3) < 1e-9


def test_build_aligned_matrix_shape():
    runs = [_make_run([0, 0, 1, 1], {}, 2), _make_run([1, 1, 0, 0], {}, 2)]
    mat, names = build_aligned_matrix(runs)
    assert mat.shape == (4, 2)
    # after alignment the two label-permuted specs agree
    assert _s3_from_matrix(mat).min() == 1.0


def test_confound_check_returns_correlation():
    rng = np.random.RandomState(0)
    labels = [0, 0, 1, 1, 2, 2]
    counts = rng.randint(1000, 5000, 6)
    runs = [_make_run(labels, {}, 3, counts=counts) for _ in range(3)]
    per_spot = pd.DataFrame({"spot_id": runs[0]["adata"].obs_names,
                             "s3_score": [1.0, 1.0, 0.6, 0.6, 0.3, 0.3]})
    res = confound_check_umi(runs, per_spot)
    assert "pearson_r" in res and -1 <= res["pearson_r"] <= 1


def test_convergence_returns_increasing_sizes():
    runs = [_make_run([0, 0, 1, 1], {}, 2) for _ in range(6)]
    df = multiverse_size_convergence(runs, sizes=[2, 4, 6], n_draws=5)
    assert list(df["n_specs"]) == [2, 4, 6]
    assert df["mean_s3"].between(0, 1).all()


def test_variance_decomposition_runs():
    runs = []
    for res in [0.3, 0.5, 0.8, 1.2] * 3:
        spec = {
            "normalization": "log1p",
            "feature_selection": {"method": "seurat", "n_genes": 3000},
            "clustering": {"method": "leiden", "resolution": res},
            "neighborhood_graph": {"type": "knn", "k_or_radius": 15},
            "dimensionality_reduction": {"n_components": 30},
        }
        # domain count driven by resolution
        runs.append(_make_run([0, 1], spec, int(5 + 10 * res)))
    df = variance_decomposition(runs)
    assert "variance_fraction" in df.columns
    assert abs(df["variance_fraction"].sum() - 1.0) < 1e-6
