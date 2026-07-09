"""Tests for src/consensus.py on a small synthetic multiverse (no network)."""

import anndata as ad
import numpy as np
import pandas as pd

from src.consensus import (
    co_association_matrix,
    consensus_labels,
    consensus_advantage_bootstrap_ci,
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


def _make_noisy_runs(n_specs=8, n=30, n_flip=3, seed=0):
    """
    Two true blocks of n/2; each spec independently mislabels a distinct random
    set of ``n_flip`` spots into the wrong block, so no single spec is perfect,
    but a spot is wrong in only a minority of specs - the co-association-based
    consensus should still recover the reference exactly.
    """
    rng = np.random.RandomState(seed)
    half = n // 2
    ref = ["A"] * half + ["B"] * (n - half)
    runs = []
    for s in range(n_specs):
        true = np.array(ref)
        labels = true.copy()
        flip_idx = rng.choice(n, n_flip, replace=False)
        labels[flip_idx] = np.where(labels[flip_idx] == "A", "B", "A")
        a = ad.AnnData(np.ones((n, 2)))
        a.obs_names = [f"s{i}" for i in range(n)]
        a.obs["leiden"] = pd.Categorical(labels.tolist())
        a.obs["cluster"] = pd.Categorical(ref)
        runs.append({"adata": a, "spec": {}, "name": f"spec{s}",
                     "metadata": {"n_clusters_found": 2}})
    return runs


def test_bootstrap_ci_zero_when_consensus_equals_best_single():
    # All specs agree perfectly -> consensus and best-single are the same
    # partition, so every bootstrap draw of the ARI difference is exactly 0.
    runs = _make_runs()
    res = consensus_advantage_bootstrap_ci(runs, k_values=[2], linkages=["average"], n_boot=200)
    assert res["observed_diff"] == 0.0
    assert res["ci95_low"] == 0.0 and res["ci95_high"] == 0.0


def test_bootstrap_ci_detects_real_consensus_advantage():
    # No single noisy spec is perfect, but the majority-vote consensus recovers
    # the true blocks exactly - a real, non-zero advantage the CI should catch.
    runs = _make_noisy_runs()
    res = evaluate_consensus_recovery(runs, k_values=[2], linkages=["average"])
    assert res["consensus_best_ari"] > res["single_max"]  # sanity: real gap exists

    boot = consensus_advantage_bootstrap_ci(runs, k_values=[2], linkages=["average"], n_boot=300)
    assert boot["observed_diff"] > 0
    assert boot["frac_bootstrap_positive"] > 0.5
    assert boot["ci95_high"] > 0  # CI doesn't collapse to a point at zero
