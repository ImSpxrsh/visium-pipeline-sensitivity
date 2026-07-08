"""
Entry point: orchestrate the full multiverse sweep.

Loads the spec config, runs every spec, saves per-spec outputs, records a failure
manifest, then computes cross-spec stability (S³) and writes aggregate results.

Run from the project root:  python scripts/main.py
"""

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

# Allow "python scripts/main.py" from the project root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.data_loader import load_dataset, prepare_data  # noqa: E402
from src.pipeline import run_pipeline  # noqa: E402
from src.stability import StabilityScorer  # noqa: E402
from src.utils import dict_to_hash, load_yaml, save_json, setup_logging  # noqa: E402

logger = logging.getLogger(__name__)


def load_specs_from_config(config_path: Path) -> Dict[str, Any]:
    """Load the config and return {'specs': [...], 'global': {...}}."""
    config = load_yaml(config_path)
    specs = config.get("specs", [])
    logger.info("Loaded %d specs from %s", len(specs), config_path)
    return config


def run_multiverse_sweep(
    adata,
    specs: List[Dict[str, Any]],
    output_dir: Path,
    random_state: int = 42,
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """Run the pipeline for each spec, saving outputs and collecting results."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results: List[Dict[str, Any]] = []
    failures: List[Dict[str, Any]] = []

    for i, spec in enumerate(specs):
        spec_name = spec.get("name", f"spec_{i:02d}")
        spec_hash = dict_to_hash(spec)
        spec_dir = output_dir / f"spec_{i:02d}_{spec_hash}"
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
                    "adata": adata_result,
                    "metadata": metadata,
                    "output_dir": str(spec_dir),
                }
            )
            logger.info("  saved -> %s", spec_dir)
        except Exception as e:  # log, record, and continue the sweep
            logger.error("  spec %s failed: %s", spec_name, e)
            failures.append({"spec_name": spec_name, "spec_hash": spec_hash, "error": str(e)})

    save_json({"failures": failures}, output_dir / "failure_manifest.json")
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
    # Persist the per-spec metadata as a specification-curve table.
    save_json(
        {"metadata": [r["metadata"] for r in results]},
        output_dir / "specification_curves.json",
    )
    logger.info("Aggregate stability saved to %s", output_dir)


def main() -> None:
    project_dir = Path(__file__).resolve().parent.parent
    results_dir = project_dir / "results"
    setup_logging(results_dir / "logs" / "sweep.log")

    logger.info("=" * 60)
    logger.info("MOSAIC: Multiverse Spatial Transcriptomics Analysis")
    logger.info("=" * 60)

    logger.info("--- Loading data ---")
    adata, _ = load_dataset("visium_hne", verbose=True)
    adata = prepare_data(adata, output_dir=project_dir / "data" / "processed", verbose=True)

    logger.info("--- Loading specifications ---")
    config = load_specs_from_config(project_dir / "configs" / "multiverse_spec.yaml")
    specs = config["specs"]
    random_state = config.get("global", {}).get("random_state", 42)

    logger.info("--- Running sweep ---")
    results = run_multiverse_sweep(
        adata, specs, output_dir=results_dir / "runs", random_state=random_state, verbose=True
    )

    logger.info("--- Computing stability ---")
    compute_aggregate_stability(results, results_root=results_dir, verbose=True)

    logger.info("=" * 60)
    logger.info("MOSAIC sweep complete")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
