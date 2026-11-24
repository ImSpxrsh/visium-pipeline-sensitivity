"""
Unit tests for mosaic/stability.py using small synthetic AnnData objects.

No network needed: we construct clusterings and Moran's I tables by hand and
check the S³ math, including the cross-spec label-alignment behavior.
"""

import anndata as ad
import numpy as np
import pandas as pd

from mosaic.stability import StabilityScorer, _align_labels_to_reference


def _make_adata(labels, svg_table=None, n_genes=3):
    n = len(labels)
    x = np.ones((n, n_genes))
    a = ad.AnnData(x)
    a.obs_names = [f"spot{i}" for i in range(n)]
    a.var_names = [f"gene{i}" for i in range(n_genes)]
    a.obs["leiden"] = pd.Categorical([str(v) for v in labels])
    a.obsm["spatial"] = np.random.RandomState(0).rand(n, 2)
    if svg_table is not None:
        a.uns["svg_results"] = svg_table
    return a


def test_align_labels_matches_permuted_labels():
    ref = np.array([0, 0, 1, 1, 2, 2])
    # Same partition, different label names -> alignment should recover ref.
    permuted = np.array([2, 2, 0, 0, 1, 1])
    aligned = _align_labels_to_reference(permuted, ref)
    assert np.array_equal(aligned, ref)


def test_per_spot_stability_all_agree():
    scorer = StabilityScorer(verbose=False)
    for _ in range(4):
        scorer.add_result({}, _make_adata([0, 0, 1, 1]), {})
    df = scorer.compute_per_spot_stability()
    assert (df["s3_score"] == 1.0).all()


def test_per_spot_stability_robust_to_label_permutation():
    scorer = StabilityScorer(verbose=False)
    scorer.add_result({}, _make_adata([0, 0, 1, 1]), {})
    scorer.add_result({}, _make_adata([1, 1, 0, 0]), {})  # same partition, swapped names
    df = scorer.compute_per_spot_stability()
    # After alignment the two specs agree perfectly.
    assert (df["s3_score"] == 1.0).all()


def test_per_spot_stability_partial_disagreement():
    scorer = StabilityScorer(verbose=False)
    scorer.add_result({}, _make_adata([0, 0, 1, 1]), {})
    scorer.add_result({}, _make_adata([0, 0, 1, 1]), {})
    scorer.add_result({}, _make_adata([0, 1, 1, 1]), {})  # spot1 flips
    df = scorer.compute_per_spot_stability().set_index("spot_id")
    assert df.loc["spot1", "s3_score"] == 2 / 3
    assert df.loc["spot0", "s3_score"] == 1.0


def test_per_feature_stability_uses_fdr_column():
    scorer = StabilityScorer(verbose=False)
    # gene0 significant in both specs, gene1 in one, gene2 in none.
    t1 = pd.DataFrame(
        {"I": [0.9, 0.5, 0.1], "pval_norm_fdr_bh": [0.001, 0.01, 0.9]},
        index=["gene0", "gene1", "gene2"],
    )
    t2 = pd.DataFrame(
        {"I": [0.8, 0.2, 0.1], "pval_norm_fdr_bh": [0.001, 0.5, 0.9]},
        index=["gene0", "gene1", "gene2"],
    )
    scorer.add_result({}, _make_adata([0, 1], svg_table=t1), {})
    scorer.add_result({}, _make_adata([0, 1], svg_table=t2), {})
    df = scorer.compute_per_feature_stability().set_index("gene")
    assert df.loc["gene0", "s3_score"] == 1.0
    assert df.loc["gene1", "s3_score"] == 0.5
    assert df.loc["gene2", "s3_score"] == 0.0


def test_partition_concordance_identical_and_disjoint():
    # Identical partitions -> ARI 1.0; label permutation must not change it.
    scorer = StabilityScorer(verbose=False)
    scorer.add_result({}, _make_adata([0, 0, 1, 1]), {})
    scorer.add_result({}, _make_adata([1, 1, 0, 0]), {})  # same partition, swapped names
    res = scorer.compute_partition_concordance()
    assert res["n_pairs"] == 1
    assert abs(res["mean"] - 1.0) < 1e-9


def test_concordance_with_reference():
    # A spec that perfectly matches the reference annotation -> ARI 1.0.
    scorer = StabilityScorer(verbose=False)
    a = _make_adata([0, 0, 1, 1])
    a.obs["cluster"] = pd.Categorical(["X", "X", "Y", "Y"])
    scorer.add_result({}, a, {})
    res = scorer.concordance_with_reference(cluster_key="leiden", reference_key="cluster")
    assert res["n_specs"] == 1
    assert abs(res["mean"] - 1.0) < 1e-9


def test_headline_instability_stats():
    scorer = StabilityScorer(verbose=False)
    for n in (5, 7, 9):
        scorer.add_result({}, _make_adata([0, 1]), {"n_clusters_found": n})
    stats = scorer.compute_headline_instability("n_clusters_found")
    assert stats["mean"] == 7.0
    assert stats["min"] == 5.0
    assert stats["max"] == 9.0
    assert stats["range"] == 4.0


def test_unmatched_clusters_do_not_share_labels():
    """Unmatched clusters get distinct labels within and across specs."""
    ref = np.array([0, 0, 1, 1, 1, 1])
    a = np.array([0, 0, 1, 1, 2, 3])  # clusters 2 and 3 cannot both match
    b = np.array([0, 0, 1, 1, 2, 3])
    la = _align_labels_to_reference(a, ref, spec_index=0)
    lb = _align_labels_to_reference(b, ref, spec_index=1)
    unmatched_a = set(la[la >= 2])
    unmatched_b = set(lb[lb >= 2])
    assert len(unmatched_a) == 2          # not merged into one overflow label
    assert unmatched_a.isdisjoint(unmatched_b)  # and never shared across specs


def test_unmatched_spots_do_not_inflate_s3():
    """Spots in unrelated unmatched clusters should not count as agreeing."""
    import anndata as ad
    import pandas as pd

    n = 6
    ref = ["0", "0", "1", "1", "1", "1"]
    other = ["0", "0", "1", "1", "2", "3"]
    sc = StabilityScorer(verbose=False)
    for i, labs in enumerate([ref, other, other]):
        a = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
        a.obs["leiden"] = pd.Categorical(labs)
        sc.add_result({"name": str(i)}, a, {})
    s3 = sc.compute_per_spot_stability("leiden")["s3_score"].values
    # The two unmatched spots never share a label across the 3 specs.
    assert np.all(s3[:4] == 1.0)
    assert np.all(s3[4:] < 1.0)
