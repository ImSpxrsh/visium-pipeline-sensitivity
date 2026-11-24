"""
Predicting fragility from a single pipeline (an early-warning shortcut).

Running the full multiverse to learn where a result is fragile is expensive. Here
we ask whether the fragile spots can be flagged from a *single* pipeline's local
structure alone -- before committing to the sweep. Features are strictly
multiverse-free (computed from one run + raw counts): local transcriptomic
heterogeneity, single-pipeline neighbor disagreement, sequencing depth, and genes
detected. A cross-validated classifier predicts which spots the full multiverse
will find fragile (S3 < 0.5).
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def single_pipeline_features(adata, k: int = 6) -> pd.DataFrame:
    """Multiverse-free per-spot features from one pipeline's output + raw counts."""
    from scipy.spatial import cKDTree

    coords = np.asarray(adata.obsm["spatial"])
    pca = np.asarray(adata.obsm["X_pca"])
    lei = adata.obs["leiden"].astype(str).values
    tree = cKDTree(coords)
    _, nbr = tree.query(coords, k=k + 1)
    nbr = nbr[:, 1:]
    het = np.array([np.linalg.norm(pca[i] - pca[nbr[i]], axis=1).mean()
                    for i in range(len(coords))])
    disagree = np.array([(lei[nbr[i]] != lei[i]).mean() for i in range(len(coords))])
    depth = (adata.obs["total_counts"].values.astype(float)
             if "total_counts" in adata.obs
             else np.asarray(adata.layers.get("counts", adata.X).sum(1)).ravel())
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    n_genes = np.asarray((counts > 0).sum(1)).ravel().astype(float)
    return pd.DataFrame({
        "local_heterogeneity": het,
        "neighbor_disagreement": disagree,
        "depth": depth,
        "n_genes": n_genes,
    }, index=adata.obs_names)


def predict_fragility(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    fragile_threshold: float = 0.5,
    n_splits: int = 5,
    seed: int = 42,
) -> Dict[str, Any]:
    """Cross-validated AUC for predicting fragile spots from a single pipeline."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    a = runs[0]["adata"]
    feats = single_pipeline_features(a)
    s3 = per_spot.set_index("spot_id").reindex(a.obs_names)["s3_score"].values
    fragile = (s3 < fragile_threshold).astype(int)

    clf = LogisticRegression(max_iter=1000)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    proba = cross_val_predict(clf, feats.values, fragile, cv=cv,
                              method="predict_proba")[:, 1]
    auc = float(roc_auc_score(fragile, proba))
    per_feat = {c: float(roc_auc_score(fragile, feats[c].values)) for c in feats.columns}
    # orient single-feature AUCs so >0.5 means "predicts fragile"
    per_feat = {c: (v if v >= 0.5 else 1 - v) for c, v in per_feat.items()}
    result = {
        "auc_combined": auc,
        "per_feature_auc": per_feat,
        "n_fragile": int(fragile.sum()),
        "n_spots": int(len(fragile)),
        "proba": proba,
        "s3": s3,
        "adata": a,
    }
    logger.info("Fragility predictable from ONE pipeline: AUC=%.3f (top feature %s=%.3f)",
                auc, max(per_feat, key=per_feat.get), max(per_feat.values()))
    return result


def spatial_block_groups(coords: np.ndarray, n_blocks: int = 10, seed: int = 42) -> np.ndarray:
    """
    Assign each spot to one of ``n_blocks`` contiguous spatial blocks via k-means
    on its (x, y) coordinates.

    Neighboring spots in Visium data are highly correlated (shared local
    heterogeneity, shared neighbor-disagreement structure), so a random
    spot-level split lets the classifier see training spots a few microns from
    its test spots -- optimistic leakage that ordinary k-fold does not control
    for. Grouping spots into spatial blocks and holding out whole blocks at a
    time removes that leakage: no block straddles both the train and test set.
    """
    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=n_blocks, random_state=seed, n_init=10)
    return km.fit_predict(coords)


def predict_fragility_buffered_spatial_cv(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    fragile_threshold: float = 0.5,
    n_blocks: int = 10,
    buffer_multiple: float = 3.0,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Buffered counterpart to ``predict_fragility_spatial_cv``.

    K-means spatial blocks remove most spatial leakage, but spots right at a
    block boundary can still be near neighbors of a held-out spot in an
    adjacent block, so "no training spot is a near neighbor of a held-out
    spot" is too strong a claim for plain spatial-block CV. This adds a
    buffer: for each fold, any training spot within ``buffer_multiple`` times
    the median nearest-neighbor spacing of a held-out spot is dropped from
    that fold's training set (not evaluated -- excluded entirely), so the
    train/test gap is at least that physical distance everywhere.
    """
    from scipy.spatial import cKDTree
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    a = runs[0]["adata"]
    feats = single_pipeline_features(a)
    s3 = per_spot.set_index("spot_id").reindex(a.obs_names)["s3_score"].values
    fragile = (s3 < fragile_threshold).astype(int)
    coords = np.asarray(a.obsm["spatial"])
    groups = spatial_block_groups(coords, n_blocks=n_blocks, seed=seed)

    tree = cKDTree(coords)
    nn_dist, _ = tree.query(coords, k=2)  # [:,0] is self (dist 0)
    buffer_radius = buffer_multiple * float(np.median(nn_dist[:, 1]))

    X = feats.values
    n = len(fragile)
    proba = np.full(n, np.nan)
    n_dropped_total = 0
    for block in np.unique(groups):
        test_mask = groups == block
        test_idx = np.where(test_mask)[0]
        candidate_train_idx = np.where(~test_mask)[0]

        # Drop any candidate-training spot within buffer_radius of a held-out spot.
        test_tree = cKDTree(coords[test_idx])
        dist_to_test, _ = test_tree.query(coords[candidate_train_idx], k=1)
        keep = dist_to_test > buffer_radius
        train_idx = candidate_train_idx[keep]
        n_dropped_total += int((~keep).sum())

        if len(np.unique(fragile[train_idx])) < 2:
            continue  # degenerate fold; leave this block's proba as NaN

        clf = LogisticRegression(max_iter=1000)
        clf.fit(X[train_idx], fragile[train_idx])
        proba[test_idx] = clf.predict_proba(X[test_idx])[:, 1]

    valid = ~np.isnan(proba)
    auc = float(roc_auc_score(fragile[valid], proba[valid]))
    result = {
        "auc_combined_buffered_spatial_cv": auc,
        "n_blocks": n_blocks,
        "buffer_radius": buffer_radius,
        "buffer_multiple": buffer_multiple,
        "mean_train_spots_dropped_per_fold": n_dropped_total / n_blocks,
        "n_fragile": int(fragile.sum()),
        "n_spots": int(len(fragile)),
        "n_spots_scored": int(valid.sum()),
        "proba": proba,
        "s3": s3,
        "valid_mask": valid,
        "adata": a,
    }
    logger.info(
        "Fragility predictable from ONE pipeline (buffered spatial-block CV): "
        "AUC=%.3f (%d blocks, buffer radius=%.1f, avg %.0f train spots dropped/fold)",
        auc, n_blocks, buffer_radius, n_dropped_total / n_blocks,
    )
    return result


def predict_fragility_spatial_cv(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    fragile_threshold: float = 0.5,
    n_blocks: int = 10,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Spatially blocked counterpart to ``predict_fragility``.

    Uses the same features, labels, and classifier as ``predict_fragility``,
    but replaces the random ``StratifiedKFold`` with group k-fold over
    spatial blocks (``spatial_block_groups``), so every fold's held-out spots
    are spatially disjoint from its training spots. Comparing this AUC to the
    one from ``predict_fragility`` quantifies how much of the original 0.81
    was genuine single-pipeline signal versus spatial-adjacency leakage.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold, cross_val_predict

    a = runs[0]["adata"]
    feats = single_pipeline_features(a)
    s3 = per_spot.set_index("spot_id").reindex(a.obs_names)["s3_score"].values
    fragile = (s3 < fragile_threshold).astype(int)
    coords = np.asarray(a.obsm["spatial"])

    groups = spatial_block_groups(coords, n_blocks=n_blocks, seed=seed)

    clf = LogisticRegression(max_iter=1000)
    cv = GroupKFold(n_splits=n_blocks)
    proba = cross_val_predict(clf, feats.values, fragile, cv=cv, groups=groups,
                              method="predict_proba")[:, 1]
    auc = float(roc_auc_score(fragile, proba))
    result = {
        "auc_combined_spatial_cv": auc,
        "n_blocks": n_blocks,
        "n_fragile": int(fragile.sum()),
        "n_spots": int(len(fragile)),
        "proba": proba,
        "s3": s3,
        "groups": groups,
        "adata": a,
    }
    logger.info("Fragility predictable from ONE pipeline (spatial-block CV): AUC=%.3f "
                "(%d spatial blocks)", auc, n_blocks)
    return result


def fragility_prediction_significance(
    fragile: np.ndarray,
    proba: np.ndarray,
    n_boot: int = 1000,
    n_perm: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Bootstrap CI and permutation null for the fragility-prediction AUC.

    ``predict_fragility`` already returns a genuine out-of-fold AUC (no
    leakage, via ``cross_val_predict``) but no uncertainty bound or null -
    CLAUDE.md requires every positive claim ship with both. Bootstraps spots
    for the CI (does not refit the classifier - ``proba`` is already valid
    held-out predictions, so only the *evaluation* is resampled) and permutes
    the fragile/not-fragile labels with ``proba`` held fixed for the null,
    matching what AUC would look like if the features carried no real signal.
    """
    from sklearn.metrics import roc_auc_score

    rng = np.random.RandomState(seed)
    n = len(fragile)
    observed = float(roc_auc_score(fragile, proba))

    boot = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(np.unique(fragile[idx])) < 2:
            continue  # degenerate resample (no positives or no negatives); skip
        boot.append(roc_auc_score(fragile[idx], proba[idx]))
    boot = np.array(boot)

    perm = np.array([
        roc_auc_score(rng.permutation(fragile), proba) for _ in range(n_perm)
    ])
    p_value = float((np.sum(perm >= observed) + 1) / (n_perm + 1))

    result = {
        "auc_observed": observed,
        "ci95_low": float(np.percentile(boot, 2.5)),
        "ci95_high": float(np.percentile(boot, 97.5)),
        "n_boot_used": int(len(boot)),
        "permutation_p_value": p_value,
        "n_perm": n_perm,
    }
    logger.info(
        "Fragility-prediction AUC significance: observed=%.3f, 95%% CI [%.3f, %.3f], "
        "permutation p=%.4f",
        observed, result["ci95_low"], result["ci95_high"], p_value,
    )
    return result


def plot_predicted_vs_actual(result: Dict[str, Any], output_path: Optional[Path] = None) -> None:
    """Predicted P(fragile) map beside the actual multiverse S3 map."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = result["adata"]
    coords = np.asarray(a.obsm["spatial"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    sc0 = axes[0].scatter(coords[:, 0], coords[:, 1], c=result["proba"], cmap="magma", s=10)
    axes[0].set_title("A  Predicted P(S3 < 0.5), one pipeline",
                      fontsize=16, fontweight="bold")
    sc1 = axes[1].scatter(coords[:, 0], coords[:, 1], c=result["s3"], cmap="RdYlGn",
                          vmin=0, vmax=1, s=10)
    axes[1].set_title("B  S3 across 30 pipelines", fontsize=16, fontweight="bold")
    for ax, sc in zip(axes, (sc0, sc1)):
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def main() -> None:
    from mosaic.utils import save_json
    from mosaic.visualization import load_runs

    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    runs = load_runs(project / "results")
    per_spot = pd.read_csv(project / "results" / "aggregate" / "per_spot_stability.csv")
    result = predict_fragility(runs, per_spot)
    fragile = (result["s3"] < 0.5).astype(int)
    significance = fragility_prediction_significance(fragile, result["proba"])
    save_json(
        {**{k: v for k, v in result.items() if k not in ("proba", "s3", "adata")},
         "significance": significance},
        project / "results" / "aggregate" / "fragility_prediction.json",
    )
    plot_predicted_vs_actual(result, project / "results" / "figures" / "17_fragility_prediction.png")

    spatial_result = predict_fragility_spatial_cv(runs, per_spot)
    spatial_fragile = (spatial_result["s3"] < 0.5).astype(int)
    spatial_significance = fragility_prediction_significance(
        spatial_fragile, spatial_result["proba"]
    )
    save_json(
        {**{k: v for k, v in spatial_result.items()
            if k not in ("proba", "s3", "adata", "groups")},
         "significance": spatial_significance,
         "random_split_auc_for_comparison": result["auc_combined"]},
        project / "results" / "aggregate" / "fragility_prediction_spatial_cv.json",
    )

    buffered = predict_fragility_buffered_spatial_cv(runs, per_spot)
    b_fragile = (buffered["s3"] < 0.5).astype(int)
    valid = ~np.isnan(buffered["proba"])
    b_significance = fragility_prediction_significance(
        b_fragile[valid], buffered["proba"][valid]
    )
    save_json(
        {**{k: v for k, v in buffered.items()
            if k not in ("proba", "s3", "adata", "groups", "valid_mask")},
         "significance": b_significance},
        project / "results" / "aggregate" / "fragility_prediction_buffered_spatial_cv.json",
    )


if __name__ == "__main__":
    main()
