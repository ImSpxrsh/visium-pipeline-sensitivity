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
    stakes_stratified_s3,
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


def test_stakes_stratified_detects_fragile_boundary():
    # A 1-D strip: left domain 0, right domain 1; the two middle spots are the
    # boundary and are given low S3, interiors high S3.
    n = 20
    a = ad.AnnData(np.ones((n, 3)))
    a.obs_names = [f"s{i}" for i in range(n)]
    a.obsm["spatial"] = np.column_stack([np.arange(n), np.zeros(n)]).astype(float)
    domains = np.array([0] * (n // 2) + [1] * (n // 2))
    runs = [{"adata": a, "spec": {}, "metadata": {}}]
    s3 = np.full(n, 0.95)
    s3[n // 2 - 1] = s3[n // 2] = 0.4  # the boundary pair is fragile
    per_spot = pd.DataFrame({"spot_id": a.obs_names, "s3_score": s3,
                             "consensus_cluster": domains})
    res, boundary = stakes_stratified_s3(runs, per_spot, k=2, n_perm=200)
    assert res["boundary_mean_s3"] < res["interior_mean_s3"]
    assert res["gap_interior_minus_boundary"] > 0
    assert res["p_value"] < 0.1


def test_fragility_report_picks_contested_gene():
    from src.statistics import fragility_report

    # gene "flip" tested in all specs, significant in 2 of 4 (freq 0.5); gene
    # "robust" always significant (freq 1.0) -> flip is the most contested.
    def svg(flip_p):
        return pd.DataFrame(
            {"I": [0.5, 0.9], "pval_norm_fdr_bh": [flip_p, 0.001]},
            index=["flip", "robust"],
        )

    runs = [
        {"adata": ad.AnnData(np.ones((3, 2))), "spec": {}, "name": f"s{i}",
         "metadata": {"spec_name": f"s{i}"}, "svg_df": svg(p)}
        for i, p in enumerate([0.01, 0.02, 0.2, 0.3])
    ]
    for r in runs:
        r["adata"].var_names = ["flip", "robust"]
    gene, s3, df = fragility_report(runs)
    assert gene == "flip"
    assert int(df["significant"].sum()) == 2


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
