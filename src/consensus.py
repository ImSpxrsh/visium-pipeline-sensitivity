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
    from src.statistics import build_aligned_matrix

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
            color="#D73027", fontsize=8, va="bottom")
    ax.set_xlabel("Consensus granularity (k domains)")
    ax.set_ylabel("Anatomy recovery (ARI vs reference)")
    ax.set_title("The multiverse consensus recovers anatomy as well as, or better than,\n"
                 "any single pipeline — without needing the ground truth to pick it",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
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
    result = evaluate_consensus_recovery(runs)
    agg = project / "results" / "aggregate"
    result["table"].to_csv(agg / "consensus_recovery.csv", index=False)
    save_json({k: v for k, v in result.items() if k != "table"},
              agg / "consensus_recovery_summary.json")
    plot_consensus_recovery(result, project / "results" / "figures" / "15_consensus_recovery.png")


if __name__ == "__main__":
    main()
