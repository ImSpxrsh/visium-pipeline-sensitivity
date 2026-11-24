"""
The multiverse as a better estimator: consensus domain recovery.

MOSAIC's core metric (S3) diagnoses instability. This module goes further: it
asks whether the *consensus* across the whole multiverse is a better estimate of
the true tissue structure than any single defensible pipeline. Because a
researcher cannot know in advance which single pipeline best matches ground truth
(that requires the ground truth), a consensus that reliably matches or beats the
best single pipeline turns the multiverse from a cautionary diagnostic into a
usable, uncertainty-quantified output.

Approach: build a co-association matrix (fraction of specifications that place
each pair of spots in the same aligned cluster), then re-cluster it. We evaluate
the consensus against the dataset's reference anatomical annotation (ARI) and
compare it to the distribution of single-pipeline ARIs.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def co_association_matrix(runs: List[Dict[str, Any]], cluster_key: str = "leiden") -> np.ndarray:
    """Fraction of specifications placing each pair of spots in the same aligned cluster."""
    from mosaic.statistics import build_aligned_matrix

    mat, _ = build_aligned_matrix(runs, cluster_key)  # spots x specs, aligned labels
    n_specs = mat.shape[1]
    n_spots = mat.shape[0]
    coassoc = np.zeros((n_spots, n_spots), dtype=np.float32)
    for s in range(n_specs):
        col = mat[:, s]
        coassoc += (col[:, None] == col[None, :])
    coassoc /= n_specs
    return coassoc


def consensus_labels(coassoc: np.ndarray, k: int, linkage: str = "average") -> np.ndarray:
    """Re-cluster the co-association matrix into k consensus domains."""
    from sklearn.cluster import AgglomerativeClustering

    dist = 1.0 - coassoc
    np.fill_diagonal(dist, 0.0)
    return AgglomerativeClustering(
        n_clusters=k, metric="precomputed", linkage=linkage
    ).fit_predict(dist)


def _reference(runs, reference_key):
    return np.array(runs[0]["adata"].obs[reference_key].astype(str).tolist())


def evaluate_consensus_recovery(
    runs: List[Dict[str, Any]],
    reference_key: str = "cluster",
    k_values: Optional[Sequence[int]] = None,
    linkages: Sequence[str] = ("average", "complete"),
    cluster_key: str = "leiden",
) -> Dict[str, Any]:
    """
    Compare consensus anatomy recovery to the single-pipeline ARI distribution.

    Returns single-pipeline ARI stats, a per-(linkage, k) consensus ARI table, and
    a summary of how often the consensus beats the best / median single pipeline.
    """
    from sklearn.metrics import adjusted_rand_score

    ref = _reference(runs, reference_key)
    n_ref = len(set(ref))
    singles = np.array([
        adjusted_rand_score(ref, np.array(r["adata"].obs[cluster_key].astype(str).tolist()))
        for r in runs
    ])
    if k_values is None:
        k_values = list(range(max(2, n_ref - 5), n_ref + 4))

    coassoc = co_association_matrix(runs, cluster_key)
    rows = []
    for linkage in linkages:
        for k in k_values:
            lab = consensus_labels(coassoc, k, linkage)
            rows.append({
                "linkage": linkage, "k": int(k),
                "ari": float(adjusted_rand_score(ref, lab)),
            })
    table = pd.DataFrame(rows)

    best_single, median_single = float(singles.max()), float(np.median(singles))
    result = {
        "n_reference_regions": n_ref,
        "single_min": float(singles.min()),
        "single_median": median_single,
        "single_max": best_single,
        "consensus_best_ari": float(table["ari"].max()),
        "consensus_at_ref_k": float(
            table[table["k"] == n_ref]["ari"].max()
        ) if (table["k"] == n_ref).any() else None,
        "beats_all_singles_frac": float((table["ari"] > best_single).mean()),
        "beats_median_frac": float((table["ari"] > median_single).mean()),
        "table": table,
    }
    logger.info(
        "Consensus recovery: single median=%.3f max=%.3f; consensus best=%.3f (at k=%d=%.3f); "
        "beats all singles in %.0f%% of settings, beats median in %.0f%%",
        median_single, best_single, result["consensus_best_ari"], n_ref,
        result["consensus_at_ref_k"] or float("nan"),
        100 * result["beats_all_singles_frac"], 100 * result["beats_median_frac"],
    )
    return result


def consensus_advantage_bootstrap_ci(
    runs: List[Dict[str, Any]],
    reference_key: str = "cluster",
    k_values: Optional[Sequence[int]] = None,
    linkages: Sequence[str] = ("average", "complete"),
    cluster_key: str = "leiden",
    n_boot: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Bootstrap CI on the ARI margin between the best consensus and the best single
    pipeline, resampling spots with replacement.

    ``evaluate_consensus_recovery`` reports point estimates only. On mouse-brain,
    the margin (consensus_best_ari - single_max) is small (~0.007 ARI) - plausibly
    within resampling noise on ~2700 spots, plausibly not. CLAUDE.md requires every
    positive claim ship with a bootstrap CI; this fills that gap for the
    consensus-estimator result rather than reporting the point estimates alone.

    Fixes the best-single and best-consensus label assignments (as already
    identified from the full dataset) and resamples the reference/predicted label
    pairs to quantify sampling variability in the ARI difference - the standard
    paired bootstrap for comparing two fixed partitions' agreement with a
    reference, not a re-derivation of the consensus itself under resampling
    (which would require re-running agglomerative clustering per draw and answers
    a different, much more expensive question).
    """
    from sklearn.metrics import adjusted_rand_score

    ref = _reference(runs, reference_key)
    n_ref = len(set(ref))
    if k_values is None:
        k_values = list(range(max(2, n_ref - 5), n_ref + 4))

    singles = [
        (np.array(r["adata"].obs[cluster_key].astype(str).tolist()), None) for r in runs
    ]
    for i, (labels, _) in enumerate(singles):
        singles[i] = (labels, float(adjusted_rand_score(ref, labels)))
    best_single_labels, single_max = max(singles, key=lambda t: t[1])

    coassoc = co_association_matrix(runs, cluster_key)
    best_ari, best_consensus_labels = -1.0, None
    for linkage in linkages:
        for k in k_values:
            lab = consensus_labels(coassoc, k, linkage)
            ari = float(adjusted_rand_score(ref, lab))
            if ari > best_ari:
                best_ari, best_consensus_labels = ari, lab

    rng = np.random.RandomState(seed)
    n = len(ref)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        diffs[b] = (
            adjusted_rand_score(ref[idx], best_consensus_labels[idx])
            - adjusted_rand_score(ref[idx], best_single_labels[idx])
        )

    result = {
        "single_max": single_max,
        "consensus_best_ari": best_ari,
        "observed_diff": best_ari - single_max,
        "bootstrap_mean_diff": float(diffs.mean()),
        "ci95_low": float(np.percentile(diffs, 2.5)),
        "ci95_high": float(np.percentile(diffs, 97.5)),
        "frac_bootstrap_positive": float((diffs > 0).mean()),
        "n_boot": n_boot,
    }
    logger.info(
        "Consensus advantage bootstrap: observed diff=%.4f, 95%% CI [%.4f, %.4f], "
        "%.0f%% of bootstrap draws favor consensus",
        result["observed_diff"], result["ci95_low"], result["ci95_high"],
        100 * result["frac_bootstrap_positive"],
    )
    return result


def s3_predicts_consensus_error(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    reference_key: str = "cluster",
    k: Optional[int] = None,
    linkage: str = "average",
    cluster_key: str = "leiden",
) -> Dict[str, Any]:
    """
    Does per-spot S3 flag where the consensus is wrong? (uncertainty quantification)

    Builds the consensus, maps each consensus domain to its majority reference
    region, and asks whether S3 separates correctly- from incorrectly-assigned
    spots (ROC AUC). A high AUC means S3 is a usable, ground-truth-free confidence
    signal for the consensus output.
    """
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score

    ref = _reference(runs, reference_key)
    k = k or len(set(ref))
    names = runs[0]["adata"].obs_names
    s3 = per_spot.set_index("spot_id").reindex(names)["s3_score"].values

    lab = consensus_labels(co_association_matrix(runs, cluster_key), k, linkage)
    pred = np.empty(len(lab), dtype=object)
    for cl in set(lab):
        m = lab == cl
        pred[m] = pd.Series(ref[m]).mode()[0]
    correct = (pred == ref)

    mask = np.isfinite(s3)
    auc = float(roc_auc_score(correct[mask], s3[mask]))
    _, p = mannwhitneyu(s3[correct & mask], s3[~correct & mask], alternative="greater")
    result = {
        "consensus_accuracy": float(correct.mean()),
        "mean_s3_correct": float(s3[correct & mask].mean()),
        "mean_s3_incorrect": float(s3[~correct & mask].mean()),
        "auc_s3_predicts_correct": auc,
        "mannwhitney_p": float(p),
    }
    logger.info("S3 predicts consensus correctness: AUC=%.3f (correct S3=%.3f vs %.3f, p=%.1e)",
                auc, result["mean_s3_correct"], result["mean_s3_incorrect"], p)
    return result


def plot_consensus_recovery(result: Dict[str, Any], output_path: Optional[Path] = None) -> None:
    """Consensus ARI vs k against the single-pipeline distribution."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    table = result["table"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axhspan(result["single_min"], result["single_max"], color="#999999", alpha=0.2,
               label="single-pipeline range")
    ax.axhline(result["single_median"], color="#555555", ls="--", label="median single pipeline")
    ax.axhline(result["single_max"], color="#333333", ls=":", label="best single pipeline")
    colors = {"average": "#1A9850", "complete": "#4575B4"}
    for linkage, sub in table.groupby("linkage"):
        ax.plot(sub["k"], sub["ari"], marker="o", color=colors.get(linkage, "#E6550D"),
                label=f"consensus ({linkage} linkage)")
    ax.axvline(result["n_reference_regions"], color="#D73027", ls="-", lw=1, alpha=0.6)
    ax.text(result["n_reference_regions"], ax.get_ylim()[0], " reference\n granularity",
            color="#D73027", fontsize=12, va="bottom")
    ax.set_xlabel("Consensus granularity (k domains)")
    ax.set_ylabel("Anatomy recovery (ARI vs reference)")
    ax.set_title("Consensus vs single pipelines (ARI with reference)",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=12, loc="lower right")
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
    result = evaluate_consensus_recovery(runs)
    agg = project / "results" / "aggregate"
    result["table"].to_csv(agg / "consensus_recovery.csv", index=False)
    per_spot = pd.read_csv(agg / "per_spot_stability.csv")
    result["uncertainty"] = s3_predicts_consensus_error(runs, per_spot)
    result["advantage_bootstrap_ci"] = consensus_advantage_bootstrap_ci(runs)
    save_json({k: v for k, v in result.items() if k != "table"},
              agg / "consensus_recovery_summary.json")
    plot_consensus_recovery(result, project / "results" / "figures" / "15_consensus_recovery.png")


if __name__ == "__main__":
    main()
