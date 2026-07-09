"""
Recompute aggregate stability (S3, ARI, headline) from saved per-spec results.

Reads results/runs/spec_*/ (written by scripts/main.py) and rewrites
results/aggregate/*. Useful to refresh aggregates after changing the stability
code without re-running the whole sweep.

Run from the project root:  PYTHONPATH=. python scripts/aggregate.py
"""

import logging
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.stability import StabilityScorer  # noqa: E402
from mosaic.utils import load_json, save_json  # noqa: E402

logger = logging.getLogger(__name__)


def load_spec_results(runs_dir: Path):
    """Yield (spec, adata, metadata) for each saved spec directory."""
    for spec_dir in sorted(runs_dir.glob("spec_*")):
        adata_path = spec_dir / "adata_clustered.h5ad"
        meta_path = spec_dir / "metadata.json"
        if not adata_path.exists() or not meta_path.exists():
            logger.warning("Skipping incomplete spec dir %s", spec_dir.name)
            continue
        adata = ad.read_h5ad(adata_path)
        # Prefer the saved SVG CSV; fall back to whatever is in uns.
        svg_csv = spec_dir / "svg_results.csv"
        if svg_csv.exists():
            adata.uns["svg_results"] = pd.read_csv(svg_csv, index_col=0)
        metadata = load_json(meta_path)
        yield metadata.get("spec", {}), adata, metadata


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    runs_dir = project / "results" / "runs"
    out_dir = project / "results" / "aggregate"
    out_dir.mkdir(parents=True, exist_ok=True)

    scorer = StabilityScorer(verbose=True)
    metadatas = []
    for spec, adata, metadata in load_spec_results(runs_dir):
        scorer.add_result(spec, adata, metadata)
        metadatas.append(metadata)

    if not scorer.results:
        logger.error("No completed spec results found in %s", runs_dir)
        sys.exit(1)

    logger.info("Aggregating %d specs", len(scorer.results))
    scorer.compute_per_spot_stability("leiden").to_csv(
        out_dir / "per_spot_stability.csv", index=False
    )
    scorer.compute_per_feature_stability().to_csv(
        out_dir / "per_feature_stability.csv", index=False
    )
    save_json(
        scorer.compute_headline_instability("n_clusters_found"),
        out_dir / "headline_instability_n_clusters.json",
    )
    save_json(
        scorer.compute_partition_concordance("leiden"),
        out_dir / "partition_concordance_ari.json",
    )
    save_json(
        scorer.concordance_with_reference("leiden", "cluster"),
        out_dir / "reference_concordance_ari.json",
    )
    save_json({"metadata": metadatas}, out_dir / "specification_curves.json")

    # Supplementary Table 1: the specification space with outcomes.
    rows = []
    for m in metadatas:
        s = m.get("spec", {})
        rows.append({
            "spec_name": m.get("spec_name", ""),
            "normalization": s.get("normalization", ""),
            "hvg_method": s.get("feature_selection", {}).get("method", ""),
            "hvg_n": s.get("feature_selection", {}).get("n_genes", ""),
            "pca_components": s.get("dimensionality_reduction", {}).get("n_components", ""),
            "graph_type": s.get("neighborhood_graph", {}).get("type", ""),
            "graph_k": s.get("neighborhood_graph", {}).get("k_or_radius", ""),
            "clustering": s.get("clustering", {}).get("method", ""),
            "resolution": s.get("clustering", {}).get("resolution", ""),
            "n_domains": m.get("n_clusters_found", ""),
            "runtime_s": round(m.get("runtime_seconds", 0) or 0, 1),
        })
    pd.DataFrame(rows).to_csv(out_dir / "supp_table1_spec_space.csv", index=False)
    logger.info("Aggregates + supplementary table written to %s", out_dir)


if __name__ == "__main__":
    main()
