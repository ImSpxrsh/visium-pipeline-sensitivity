"""
Robustness of per-spot S3 to the Hungarian-alignment reference spec.

compute_per_spot_stability aligns every spec's cluster labels to one reference
spec via Hungarian (maximum-overlap) assignment before computing S3. Because
cluster *labels* are arbitrary, a different reference should not materially
change which spots are called fragile. This module checks that directly: it
recomputes per-spot S3 once per candidate reference and reports how much the
resulting S3 vectors agree.
"""

import logging
from typing import Any, Dict, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def reference_robustness(runs: List[Dict[str, Any]], cluster_key: str = "leiden") -> Dict[str, Any]:
    """
    Recompute per-spot S3 once per spec used as the alignment reference.

    Returns pairwise Pearson/Spearman correlation statistics across the
    resulting S3 vectors, plus the per-spot standard deviation across
    reference choices (how much a single spot's S3 moves depending on which
    spec anchors the alignment).
    """
    from scipy.stats import pearsonr, spearmanr

    from mosaic.stability import StabilityScorer

    n = len(runs)
    s3_by_ref = {}
    for i in range(n):
        scorer = StabilityScorer(verbose=False)
        for r in runs:
            scorer.add_result(r["spec"], r["adata"], r["metadata"])
        df = scorer.compute_per_spot_stability(cluster_key=cluster_key, reference_idx=i)
        s3_by_ref[i] = df.set_index("spot_id")["s3_score"]
        logger.info("Reference spec %d/%d (%s): mean S3=%.4f", i + 1, n,
                    runs[i]["name"], df["s3_score"].mean())

    mat = pd.DataFrame(s3_by_ref)  # spots x n_references
    spot_std = mat.std(axis=1)
    spot_mean = mat.mean(axis=1)

    pearson_rs, spearman_rs = [], []
    for i in range(n):
        for j in range(i + 1, n):
            pearson_rs.append(pearsonr(mat[i], mat[j])[0])
            spearman_rs.append(spearmanr(mat[i], mat[j])[0])

    result = {
        "n_references": n,
        "n_spots": int(mat.shape[0]),
        "mean_pairwise_pearson_r": float(np.mean(pearson_rs)),
        "min_pairwise_pearson_r": float(np.min(pearson_rs)),
        "mean_pairwise_spearman_r": float(np.mean(spearman_rs)),
        "min_pairwise_spearman_r": float(np.min(spearman_rs)),
        "mean_per_spot_std_across_refs": float(spot_std.mean()),
        "median_per_spot_std_across_refs": float(spot_std.median()),
        "p95_per_spot_std_across_refs": float(spot_std.quantile(0.95)),
        "mean_s3_range_across_refs": float((mat.max(axis=1) - mat.min(axis=1)).mean()),
        "s3_across_refs_mean": float(spot_mean.mean()),
        "s3_across_refs_std_of_means": float(
            mat.mean(axis=0).std()
        ),  # spread of each reference's own overall mean S3
    }
    logger.info(
        "Reference robustness: mean pairwise Pearson r=%.3f (min %.3f), "
        "mean per-spot std across %d references=%.4f",
        result["mean_pairwise_pearson_r"], result["min_pairwise_pearson_r"],
        n, result["mean_per_spot_std_across_refs"],
    )
    return result


def main() -> None:
    import json
    from pathlib import Path

    from mosaic.utils import save_json
    from mosaic.visualization import load_runs

    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    runs = load_runs(project / "results")
    result = reference_robustness(runs)
    save_json(result, project / "results" / "aggregate" / "reference_robustness.json")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
