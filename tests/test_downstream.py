"""Tests for mosaic/downstream.py on a small synthetic sweep (no network)."""

import anndata as ad
import numpy as np
import pandas as pd

from mosaic.downstream import _brier_ece, downstream_region_classification, summarize


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


def test_brier_ece_perfect_and_worst():
    classes = np.array(["A", "B"])
    y = np.array(["A", "B", "A", "B"])
    perfect = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
    brier, ece = _brier_ece(y, perfect, classes)
    assert brier < 1e-9 and ece < 1e-9
    # confidently wrong -> Brier = 2.0 per sample
    wrong = np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
    brier_w, ece_w = _brier_ece(y, wrong, classes)
    assert abs(brier_w - 2.0) < 1e-9 and abs(ece_w - 1.0) < 1e-9


def test_summarize_keys():
    runs = [_make_run(s) for s in range(3)]
    s = summarize(downstream_region_classification(runs, n_splits=3))
    assert {"n_specs", "accuracy_min", "accuracy_max", "accuracy_range_pp"} <= set(s)
