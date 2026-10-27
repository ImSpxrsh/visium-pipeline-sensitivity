"""
Ground-truth-free selection of the consensus (linkage, k).

evaluate_consensus_recovery (mosaic/consensus.py) sweeps (linkage, k) and reports
the ARI against the reference annotation for every setting -- useful for showing
the ceiling, but it does not show how a researcher without ground truth would
pick a single setting. This module selects (linkage, k) using only the
co-association matrix itself (mean internal silhouette score, no reference
labels) and reports how close that blind pick comes to the oracle-best ARI.
"""

import logging
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def select_consensus_ground_truth_free(
    runs: List[Dict[str, Any]],
    k_values: Optional[Sequence[int]] = None,
    linkages: Sequence[str] = ("average", "complete"),
    cluster_key: str = "leiden",
) -> Dict[str, Any]:
    """
    Pick (linkage, k) by mean silhouette on the co-association distance alone
    (no reference labels), then report that pick's actual ARI against the
    reference next to the oracle-best ARI from the same (linkage, k) grid.
    """
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    from mosaic.consensus import _reference, co_association_matrix, consensus_labels

    ref = _reference(runs, "cluster")
    n_ref = len(set(ref))
    if k_values is None:
        k_values = list(range(max(2, n_ref - 5), n_ref + 4))

    coassoc = co_association_matrix(runs, cluster_key)
    dist = 1.0 - coassoc
    np.fill_diagonal(dist, 0.0)

    rows = []
    for linkage in linkages:
        for k in k_values:
            lab = consensus_labels(coassoc, k, linkage)
            if len(set(lab)) < 2 or len(set(lab)) >= len(lab):
                continue
            sil = float(silhouette_score(dist, lab, metric="precomputed"))
            ari = float(adjusted_rand_score(ref, lab))
            rows.append({"linkage": linkage, "k": int(k), "silhouette": sil, "ari": ari})
    table = pd.DataFrame(rows)

    blind_pick = table.loc[table["silhouette"].idxmax()]
    oracle_pick = table.loc[table["ari"].idxmax()]

    result = {
        "blind_linkage": blind_pick["linkage"],
        "blind_k": int(blind_pick["k"]),
        "blind_silhouette": float(blind_pick["silhouette"]),
        "blind_ari": float(blind_pick["ari"]),
        "oracle_linkage": oracle_pick["linkage"],
        "oracle_k": int(oracle_pick["k"]),
        "oracle_ari": float(oracle_pick["ari"]),
        "ari_gap_blind_vs_oracle": float(oracle_pick["ari"] - blind_pick["ari"]),
        "silhouette_ari_spearman": float(
            table["silhouette"].corr(table["ari"], method="spearman")
        ),
        "n_settings_evaluated": int(len(table)),
        "table": table,
    }
    logger.info(
        "Ground-truth-free consensus pick: %s/k=%d (silhouette=%.3f) -> ARI=%.3f; "
        "oracle best %s/k=%d -> ARI=%.3f (gap=%.3f); silhouette-ARI Spearman rho=%.3f",
        result["blind_linkage"], result["blind_k"], result["blind_silhouette"],
        result["blind_ari"], result["oracle_linkage"], result["oracle_k"],
        result["oracle_ari"], result["ari_gap_blind_vs_oracle"],
        result["silhouette_ari_spearman"],
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
    result = select_consensus_ground_truth_free(runs)
    save_json(
        {k: v for k, v in result.items() if k != "table"},
        project / "results" / "aggregate" / "consensus_selection_blind.json",
    )
    result["table"].to_csv(
        project / "results" / "aggregate" / "consensus_selection_blind_table.csv", index=False
    )
    print(json.dumps({k: v for k, v in result.items() if k != "table"}, indent=2))


if __name__ == "__main__":
    main()
