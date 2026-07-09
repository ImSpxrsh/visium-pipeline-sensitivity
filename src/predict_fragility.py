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
    axes[0].set_title(f"A  Predicted P(fragile) from ONE pipeline\n(AUC {result['auc_combined']:.2f})",
                      fontsize=12, fontweight="bold")
    sc1 = axes[1].scatter(coords[:, 0], coords[:, 1], c=result["s3"], cmap="RdYlGn",
                          vmin=0, vmax=1, s=10)
    axes[1].set_title("B  Actual multiverse S3 (30 pipelines)", fontsize=12, fontweight="bold")
    for ax, sc in zip(axes, (sc0, sc1)):
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Fragility is predictable before running the multiverse",
                 fontsize=13, fontweight="bold")
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def main() -> None:
    from src.utils import save_json
    from src.visualization import load_runs

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


if __name__ == "__main__":
    main()
