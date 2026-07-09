"""Tests for src/consensus.py on a small synthetic multiverse (no network)."""

import anndata as ad
import numpy as np
import pandas as pd

from src.consensus import (
    co_association_matrix,
    consensus_labels,
    evaluate_consensus_recovery,
)


def _make_runs(n_specs=6, n=12):
    """All specs agree on two blocks (0..n/2, n/2..n); reference matches the blocks."""
    runs = []
    half = n // 2
    ref = ["A"] * half + ["B"] * (n - half)
    for s in range(n_specs):
        a = ad.AnnData(np.ones((n, 2)))
        a.obs_names = [f"s{i}" for i in range(n)]
        # arbitrary label names per spec (tests alignment), same partition
        labels = ([str(s)] * half) + ([str(s + 100)] * (n - half))
        a.obs["leiden"] = pd.Categorical(labels)
        a.obs["cluster"] = pd.Categorical(ref)
        runs.append({"adata": a, "spec": {}, "name": f"spec{s}",
                     "metadata": {"n_clusters_found": 2}})
    return runs


def test_co_association_is_block_when_all_agree():
    runs = _make_runs()
    c = co_association_matrix(runs)
    n = c.shape[0]
    half = n // 2
    # within-block pairs co-associate fully; cross-block pairs never
    assert c[:half, :half].min() == 1.0
    assert c[:half, half:].max() == 0.0


def test_consensus_recovers_blocks_and_beats_nothing_when_all_agree():
    runs = _make_runs()
    c = co_association_matrix(runs)
    lab = consensus_labels(c, k=2, linkage="average")
    # two recovered groups match the two reference blocks
    assert len(set(lab)) == 2
    res = evaluate_consensus_recovery(runs, k_values=[2], linkages=["average"])
    assert res["consensus_best_ari"] == 1.0  # perfect recovery
    assert res["single_max"] == 1.0
