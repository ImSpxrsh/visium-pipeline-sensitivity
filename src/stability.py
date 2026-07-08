"""
Compute the Spatial Stability Score (S³) across the multiverse.

Three views of stability:
    * per-spot   -> how often a spot lands in its consensus domain
    * per-feature-> how often a gene is called a significant SVG
    * headline   -> spread of a summary metric (e.g. cluster count) across specs

Cross-spec caveat: cluster *labels* are arbitrary (cluster "3" in one spec is
unrelated to cluster "3" in another). Before comparing labels we align every
spec's labels to a reference spec by maximum-overlap (Hungarian) matching. This
makes per-spot S³ meaningful rather than an artifact of label numbering.
"""

import logging
from typing import Any, Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _align_labels_to_reference(
    labels: np.ndarray, ref_labels: np.ndarray
) -> np.ndarray:
    """
    Relabel ``labels`` so its clusters best match ``ref_labels`` by overlap.

    Uses the Hungarian algorithm on the negative contingency table. Returns an
    integer array aligned to the reference label space.
    """
    from scipy.optimize import linear_sum_assignment

    ref_codes = pd.Categorical(ref_labels).codes
    cur = pd.Categorical(labels)
    cur_codes = cur.codes

    n_ref = ref_codes.max() + 1
    n_cur = cur_codes.max() + 1
    contingency = np.zeros((n_cur, n_ref), dtype=int)
    for c, r in zip(cur_codes, ref_codes):
        contingency[c, r] += 1

    # Maximize overlap -> minimize negative overlap.
    row_ind, col_ind = linear_sum_assignment(-contingency)
    mapping = {int(r): int(c) for r, c in zip(row_ind, col_ind)}

    # Clusters with no assigned reference (n_cur > n_ref) keep a fresh label.
    next_label = n_ref
    aligned = np.empty_like(cur_codes)
    for i, code in enumerate(cur_codes):
        if code in mapping:
            aligned[i] = mapping[code]
        else:
            aligned[i] = next_label
    return aligned


class StabilityScorer:
    """Accumulate per-spec results and compute S³ across them."""

    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.results: List[Tuple[Dict[str, Any], ad.AnnData, Dict[str, Any]]] = []

    def add_result(
        self, spec: Dict[str, Any], adata: ad.AnnData, metadata: Dict[str, Any]
    ) -> None:
        """Register one pipeline result (spec, clustered AnnData, metadata)."""
        self.results.append((spec, adata, metadata))

    def compute_per_spot_stability(self, cluster_key: str = "leiden") -> pd.DataFrame:
        """
        Per-spot S³: fraction of specs assigning each spot to its consensus domain.

        Labels are aligned to the first spec before comparison. Returns a DataFrame
        with columns [spot_id, s3_score, consensus_cluster].
        """
        if not self.results:
            raise ValueError("No results to compute stability on")

        # Collect (aligned) label vectors from specs that share the reference index.
        ref_adata = self.results[0][1]
        ref_labels = ref_adata.obs[cluster_key].values

        aligned_cols = []
        for _, adata, _ in self.results:
            if cluster_key not in adata.obs:
                logger.warning("Missing cluster key %r; skipping a spec", cluster_key)
                continue
            if not adata.obs_names.equals(ref_adata.obs_names):
                logger.warning("Spot index mismatch vs reference; skipping a spec")
                continue
            aligned_cols.append(
                _align_labels_to_reference(adata.obs[cluster_key].values, ref_labels)
            )

        cluster_matrix = np.array(aligned_cols).T  # (n_spots, n_specs_used)
        n_specs = cluster_matrix.shape[1]
        if self.verbose:
            logger.info("Per-spot cluster matrix: %s", cluster_matrix.shape)

        s3_scores = np.empty(cluster_matrix.shape[0])
        consensus = np.empty(cluster_matrix.shape[0], dtype=int)
        for i in range(cluster_matrix.shape[0]):
            row = cluster_matrix[i]
            values, counts = np.unique(row, return_counts=True)
            top = counts.argmax()
            consensus[i] = int(values[top])
            s3_scores[i] = counts[top] / n_specs

        df = pd.DataFrame(
            {
                "spot_id": ref_adata.obs_names,
                "s3_score": s3_scores,
                "consensus_cluster": consensus,
            }
        )
        if self.verbose:
            logger.info(
                "Per-spot S3: mean=%.3f std=%.3f over %d spots",
                df["s3_score"].mean(),
                df["s3_score"].std(),
                len(df),
            )
        return df

    def compute_per_feature_stability(
        self, svg_key: str = "svg_results", alpha: float = 0.05
    ) -> pd.DataFrame:
        """
        Per-feature S³: fraction of specs where a gene is a significant SVG.

        Significance uses the FDR-corrected p-value column from squidpy's Moran's I
        output when available, falling back to the raw p-value, then to an
        above-median Moran's I heuristic. Returns [gene, s3_score, n_specs_present].
        """
        if not self.results:
            raise ValueError("No results to compute stability on")

        n_specs = len(self.results)
        svg_calls: Dict[str, List[bool]] = {}

        for _, adata, _ in self.results:
            if svg_key not in adata.uns:
                logger.warning("Missing SVG key %r; skipping a spec", svg_key)
                continue
            svg_df = adata.uns[svg_key]
            is_svg = self._significant_svgs(svg_df, alpha)
            for gene, sig in zip(svg_df.index.tolist(), is_svg):
                svg_calls.setdefault(gene, []).append(bool(sig))

        rows = []
        for gene, calls in svg_calls.items():
            rows.append(
                {
                    "gene": gene,
                    "s3_score": float(np.sum(calls) / n_specs),
                    "n_specs_present": len(calls),
                }
            )
        df = pd.DataFrame(rows).sort_values("s3_score", ascending=False).reset_index(drop=True)

        if self.verbose and len(df):
            logger.info(
                "Per-feature S3: robust (>0.8)=%d, fragile (<0.5)=%d, of %d genes",
                int((df["s3_score"] > 0.8).sum()),
                int((df["s3_score"] < 0.5).sum()),
                len(df),
            )
        return df

    @staticmethod
    def _significant_svgs(svg_df: pd.DataFrame, alpha: float) -> np.ndarray:
        """Pick the best available significance signal from a Moran's I table."""
        for col in ("pval_norm_fdr_bh", "pval_sim_fdr_bh", "pval_norm", "pval_sim", "pval"):
            if col in svg_df.columns:
                return (svg_df[col] < alpha).values
        # No p-value column: fall back to a top-half Moran's I heuristic.
        return (svg_df["I"] > svg_df["I"].median()).values

    def compute_headline_instability(self, metric: str = "n_clusters_found") -> Dict[str, Any]:
        """Summarize the spread of a scalar metadata metric across specs."""
        values = [m[metric] for _, _, m in self.results if metric in m]
        if not values:
            logger.warning("No values found for metric %r", metric)
            return {}
        arr = np.asarray(values, dtype=float)
        result = {
            "metric": metric,
            "n_specs": int(arr.size),
            "mean": float(arr.mean()),
            "std": float(arr.std()),
            "min": float(arr.min()),
            "max": float(arr.max()),
            "range": float(arr.max() - arr.min()),
        }
        if self.verbose:
            logger.info(
                "Headline %s: %.2f +/- %.2f, range %.0f-%.0f",
                metric, result["mean"], result["std"], result["min"], result["max"],
            )
        return result

    def compute_partition_concordance(self, cluster_key: str = "leiden") -> Dict[str, Any]:
        """
        Label-invariant clustering agreement across specs via Adjusted Rand Index.

        For every pair of specifications we compute the ARI between their spot
        partitions (ARI needs no label alignment: 1.0 = identical partitions,
        ~0.0 = chance agreement). Returns the mean/min/max over all pairs, a
        standard and reference-free complement to the aligned per-spot S3.
        """
        from itertools import combinations

        from sklearn.metrics import adjusted_rand_score

        ref = self.results[0][1]
        labels = [
            adata.obs[cluster_key].values
            for _, adata, _ in self.results
            if cluster_key in adata.obs and adata.obs_names.equals(ref.obs_names)
        ]
        if len(labels) < 2:
            return {}

        aris = [adjusted_rand_score(a, b) for a, b in combinations(labels, 2)]
        aris = np.asarray(aris)
        result = {
            "metric": "adjusted_rand_index",
            "n_pairs": int(aris.size),
            "mean": float(aris.mean()),
            "std": float(aris.std()),
            "min": float(aris.min()),
            "max": float(aris.max()),
        }
        if self.verbose:
            logger.info(
                "Partition concordance (ARI): mean=%.3f over %d pairs (min=%.3f, max=%.3f)",
                result["mean"], result["n_pairs"], result["min"], result["max"],
            )
        return result

    def concordance_with_reference(
        self, cluster_key: str = "leiden", reference_key: str = "cluster"
    ) -> Dict[str, Any]:
        """
        Agreement (ARI) between each spec's partition and a reference annotation.

        Uses the anatomical annotation shipped with the dataset (e.g. brain
        regions in ``adata.obs['cluster']``) as ground truth. Returns the mean/
        min/max ARI over specs, quantifying how well the multiverse recovers
        known structure.
        """
        from sklearn.metrics import adjusted_rand_score

        aris = []
        for _, adata, _ in self.results:
            if cluster_key in adata.obs and reference_key in adata.obs:
                aris.append(
                    adjusted_rand_score(adata.obs[reference_key], adata.obs[cluster_key])
                )
        if not aris:
            return {}
        aris = np.asarray(aris)
        result = {
            "metric": "ari_vs_reference",
            "reference_key": reference_key,
            "n_specs": int(aris.size),
            "mean": float(aris.mean()),
            "std": float(aris.std()),
            "min": float(aris.min()),
            "max": float(aris.max()),
        }
        if self.verbose:
            logger.info(
                "Concordance with reference (ARI): mean=%.3f over %d specs (min=%.3f, max=%.3f)",
                result["mean"], result["n_specs"], result["min"], result["max"],
            )
        return result

    def summary_report(self) -> Dict[str, Any]:
        """Bundle the stability views into one dict."""
        return {
            "n_specs": len(self.results),
            "per_spot_stability": self.compute_per_spot_stability(),
            "per_feature_stability": self.compute_per_feature_stability(),
            "partition_concordance": self.compute_partition_concordance(),
            "headline_instability": {
                "n_clusters": self.compute_headline_instability("n_clusters_found"),
                "runtime": self.compute_headline_instability("runtime_seconds"),
            },
        }


def main() -> None:
    """Run two specs and compute stability between them as a smoke test."""
    from src.data_loader import load_dataset, prepare_data
    from src.pipeline import run_pipeline

    logging.basicConfig(level=logging.INFO)
    adata, _ = load_dataset("visium_hne", verbose=True)
    adata = prepare_data(adata, verbose=True)

    base = {
        "normalization": "log1p",
        "feature_selection": {"method": "seurat", "n_genes": 3000},
        "dimensionality_reduction": {"n_components": 30},
        "neighborhood_graph": {"type": "knn", "k_or_radius": 10},
        "clustering": {"method": "leiden", "resolution": 0.5},
        "spatially_variable_genes": {"method": "morans_i"},
    }
    fine = {**base, "clustering": {"method": "leiden", "resolution": 0.8}}

    scorer = StabilityScorer(verbose=True)
    for spec in (base, fine):
        r = run_pipeline(adata.copy(), spec, verbose=True)
        scorer.add_result(spec, r["adata"], r["metadata"])

    scorer.summary_report()
    logger.info("Stability computation complete")


if __name__ == "__main__":
    main()
