"""Tests for src/downstream.py on a small synthetic sweep (no network)."""

import anndata as ad
import numpy as np
import pandas as pd

from src.downstream import downstream_region_classification, summarize


def _make_run(seed, n=60, d=5):
    rng = np.random.RandomState(seed)
    # Two well-separated region clusters in PCA space -> classifier should do well.
    half = n // 2
    pca = np.vstack([rng.normal(-2, 0.5, (half, d)), rng.normal(2, 0.5, (n - half, d))])
    a = ad.AnnData(np.ones((n, 3)))
    a.obs_names = [f"s{i}" for i in range(n)]
    a.obsm["X_pca"] = pca
    a.obs["cluster"] = pd.Categorical(["A"] * half + ["B"] * (n - half))
    a.obs["leiden"] = pd.Categorical(["0"] * n)
    return {"adata": a, "spec": {"clustering": {"resolution": 0.5}},
            "metadata": {"spec_name": f"spec{seed}", "n_clusters_found": 2}}


def test_downstream_runs_and_scores_high():
    runs = [_make_run(s) for s in range(4)]
    df = downstream_region_classification(runs, n_splits=3)
    assert len(df) == 4
    assert df["accuracy"].between(0, 1).all()
    # separable clusters -> high accuracy
    assert df["accuracy"].min() > 0.8


def test_summarize_keys():
    runs = [_make_run(s) for s in range(3)]
    s = summarize(downstream_region_classification(runs, n_splits=3))
    assert {"n_specs", "accuracy_min", "accuracy_max", "accuracy_range_pp"} <= set(s)
