"""Test src/predict_fragility.py on a synthetic single pipeline (no network)."""

import anndata as ad
import numpy as np
import pandas as pd

from src.predict_fragility import predict_fragility, single_pipeline_features


def _make_run(n=120, d=5):
    rng = np.random.RandomState(0)
    # two spatial blocks on a line; the middle (boundary) spots are fragile
    coords = np.column_stack([np.arange(n), np.zeros(n)]).astype(float)
    leiden = np.array(["A"] * (n // 2) + ["B"] * (n - n // 2))
    a = ad.AnnData(rng.poisson(3, (n, d)).astype(float))
    a.obs_names = [f"s{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(d)]
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = coords
    a.obsm["X_pca"] = rng.normal(0, 1, (n, d))
    a.obs["leiden"] = pd.Categorical(leiden)
    a.obs["total_counts"] = np.asarray(a.X.sum(1)).ravel()
    return {"adata": a, "spec": {}, "name": "s0", "metadata": {}}


def test_features_shape_and_columns():
    run = _make_run()
    f = single_pipeline_features(run["adata"])
    assert set(f.columns) == {"local_heterogeneity", "neighbor_disagreement", "depth", "n_genes"}
    assert len(f) == run["adata"].n_obs


def test_predict_fragility_learns_boundary():
    run = _make_run()
    n = run["adata"].n_obs
    # fragile exactly at the two boundary spots (high neighbor disagreement)
    s3 = np.full(n, 0.95)
    s3[n // 2 - 1] = s3[n // 2] = 0.3
    # give a few more fragile spots near the boundary so CV has positives per fold
    s3[n // 2 - 3: n // 2 + 3] = 0.3
    per_spot = pd.DataFrame({"spot_id": run["adata"].obs_names, "s3_score": s3})
    res = predict_fragility([run], per_spot, n_splits=3)
    assert 0.0 <= res["auc_combined"] <= 1.0
    assert res["per_feature_auc"]["neighbor_disagreement"] >= 0.5
