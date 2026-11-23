"""
Entry point: orchestrate the full multiverse sweep.

Loads the spec config, runs every spec, saves per-spec outputs, records a failure
manifest, then computes cross-spec stability (S³) and writes aggregate results.

Run from the project root:  python scripts/main.py
"""

import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Allow "python scripts/main.py" from the project root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.data_loader import load_dataset, prepare_data  # noqa: E402
from mosaic.pipeline import run_pipeline  # noqa: E402
from mosaic.stability import StabilityScorer  # noqa: E402
from mosaic.utils import dict_to_hash, load_yaml, save_json, setup_logging  # noqa: E402

logger = logging.getLogger(__name__)


def load_specs_from_config(config_path: Path) -> Dict[str, Any]:
    """Load the config and return {'specs': [...], 'global': {...}}."""
    config = load_yaml(config_path)
    specs = config.get("specs", [])
    logger.info("Loaded %d specs from %s", len(specs), config_path)
    return config


def _slim(adata):
    """Keep only what the stability scorer reads (obs labels, obsm, the SVG table).

    Dropping X, layers, and neighbor graphs cuts memory from ~150 MB to a few MB
    per spec, so hundreds of specs fit in memory during aggregation.
    """
    import anndata as ad

    return ad.AnnData(
        obs=adata.obs.copy(),
        obsm={k: v for k, v in adata.obsm.items()},
        # Only the SVG table: uns also holds the tissue images (~50 MB per spec).
        uns={k: adata.uns[k] for k in ("svg_results",) if k in adata.uns},
    )


def run_multiverse_sweep(
    adata,
    specs: List[Dict[str, Any]],
    output_dir: Path,
    random_state: int = 42,
    verbose: bool = True,
    only_index: Optional[int] = None,
    lean: bool = False,
) -> List[Dict[str, Any]]:
    """Run the pipeline for each spec, saving outputs and collecting results.

    If ``only_index`` is given, run just that spec (used by cluster job arrays,
    one spec per job); its output directory keeps the same index and hash.
    If ``lean`` is true, the returned AnnData objects keep only obs, obsm, and
    uns (see ``_slim``); saved files are unaffected.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results: List[Dict[str, Any]] = []
    failures: List[Dict[str, Any]] = []

    for i, spec in enumerate(specs):
        if only_index is not None and i != only_index:
            continue
        spec_name = spec.get("name", f"spec_{i:02d}")
        spec_hash = dict_to_hash(spec)
        spec_dir = output_dir / f"spec_{i:02d}_{spec_hash}"

        # Incremental execution: skip a spec whose output (matched by content hash,
        # not index) already exists, so re-running only computes new specs.
        existing = next(
            (p for p in output_dir.glob(f"spec_*_{spec_hash}")
             if (p / "adata_clustered.h5ad").exists()),
            None,
        )
        if existing is not None:
            import anndata as ad

            logger.info("[%d/%d] %s already computed -> %s (skip)",
                        i + 1, len(specs), spec_name, existing.name)
            all_results.append({
                "spec": spec,
                "adata": (_slim if lean else (lambda a: a))(
                    ad.read_h5ad(existing / "adata_clustered.h5ad")
                ),
                "metadata": json.loads((existing / "metadata.json").read_text()),
                "output_dir": str(existing),
            })
            continue

        spec_dir.mkdir(parents=True, exist_ok=True)
        logger.info("[%d/%d] Running %s ...", i + 1, len(specs), spec_name)
        try:
            result = run_pipeline(adata.copy(), spec, random_state=random_state, verbose=verbose)
            adata_result = result["adata"]
            adata_result.write_h5ad(spec_dir / "adata_clustered.h5ad")

            if result["svg_results"] is not None:
                result["svg_results"].to_csv(spec_dir / "svg_results.csv")

            metadata = result["metadata"]
            metadata["spec_name"] = spec_name
            metadata["spec_hash"] = spec_hash
            save_json(metadata, spec_dir / "metadata.json")

            all_results.append(
                {
                    "spec": spec,
                    "adata": _slim(adata_result) if lean else adata_result,
                    "metadata": metadata,
                    "output_dir": str(spec_dir),
                }
            )
            logger.info("  saved -> %s", spec_dir)
        except Exception as e:  # log, record, and continue the sweep
            logger.error("  spec %s failed: %s", spec_name, e)
            failures.append({"spec_name": spec_name, "spec_hash": spec_hash, "error": str(e)})

    # Parallel array jobs each write their own manifest so they don't overwrite each other.
    manifest = ("failure_manifest.json" if only_index is None
                else f"failure_manifest_spec{only_index}.json")
    save_json({"failures": failures}, output_dir / manifest)
    logger.info(
        "Sweep complete: %d/%d specs succeeded (%d failed)",
        len(all_results), len(specs), len(failures),
    )
    return all_results


def compute_aggregate_stability(
    results: List[Dict[str, Any]], results_root: Path, verbose: bool = True
) -> None:
    """Compute S³ across specs and write aggregate CSV/JSON outputs."""
    if not results:
        logger.warning("No successful results; skipping stability computation")
        return

    output_dir = Path(results_root) / "aggregate"
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Computing stability across %d specs...", len(results))

    scorer = StabilityScorer(verbose=verbose)
    for r in results:
        scorer.add_result(r["spec"], r["adata"], r["metadata"])

    scorer.compute_per_spot_stability(cluster_key="leiden").to_csv(
        output_dir / "per_spot_stability.csv", index=False
    )
    scorer.compute_per_feature_stability().to_csv(
        output_dir / "per_feature_stability.csv", index=False
    )
    save_json(
        scorer.compute_headline_instability("n_clusters_found"),
        output_dir / "headline_instability_n_clusters.json",
    )
    save_json(
        scorer.compute_partition_concordance(cluster_key="leiden"),
        output_dir / "partition_concordance_ari.json",
    )
    save_json(
        scorer.concordance_with_reference("leiden", "cluster"),
        output_dir / "reference_concordance_ari.json",
    )
    # Persist the per-spec metadata as a specification-curve table.
    save_json(
        {"metadata": [r["metadata"] for r in results]},
        output_dir / "specification_curves.json",
    )
    logger.info("Aggregate stability saved to %s", output_dir)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the MOSAIC multiverse sweep.")
    parser.add_argument(
        "--dataset", default="visium_hne",
        help="squidpy dataset name (e.g. visium_hne, visium_fluo)",
    )
    parser.add_argument(
        "--results-root", default=None,
        help="output root (default: results/ for visium_hne, results/<dataset>/ otherwise)",
    )
    parser.add_argument(
        "--config", default=None,
        help="spec file (default: configs/multiverse_spec.yaml)",
    )
    parser.add_argument(
        "--spec-index", type=int, default=None,
        help="run only the spec at this 0-based index and skip aggregation "
             "(for cluster job arrays; rerun without it afterwards to aggregate)",
    )
    parser.add_argument(
        "--lean", action="store_true",
        help="keep only labels and SVG tables in memory while aggregating "
             "(needed for sweeps with hundreds of specs)",
    )
    args = parser.parse_args()

    project_dir = Path(__file__).resolve().parent.parent
    # Keep the primary dataset at results/; namespace any other dataset.
    if args.results_root is not None:
        results_dir = Path(args.results_root)
    elif args.dataset == "visium_hne":
        results_dir = project_dir / "results"
    else:
        results_dir = project_dir / "results" / args.dataset
    setup_logging(results_dir / "logs" / "sweep.log")

    logger.info("=" * 60)
    logger.info("MOSAIC: Multiverse Spatial Transcriptomics Analysis (%s)", args.dataset)
    logger.info("=" * 60)

    logger.info("--- Loading data ---")
    adata, _ = load_dataset(args.dataset, verbose=True)
    adata = prepare_data(adata, verbose=True)

    logger.info("--- Loading specifications ---")
    config_path = Path(args.config) if args.config else project_dir / "configs" / "multiverse_spec.yaml"
    config = load_specs_from_config(config_path)
    specs = config["specs"]
    random_state = config.get("global", {}).get("random_state", 42)

    logger.info("--- Running sweep ---")
    results = run_multiverse_sweep(
        adata, specs, output_dir=results_dir / "runs", random_state=random_state,
        verbose=True, only_index=args.spec_index, lean=args.lean,
    )
    if args.spec_index is not None:
        logger.info("Ran spec %d only; skipping aggregation", args.spec_index)
        return

    logger.info("--- Computing stability ---")
    compute_aggregate_stability(results, results_root=results_dir, verbose=True)

    logger.info("=" * 60)
    logger.info("MOSAIC sweep complete")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
