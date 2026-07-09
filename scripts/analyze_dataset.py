"""
Run the MOSAIC analysis suite on one dataset's completed sweep.

Produces the statistics (variance decomposition [H2], stakes-stratified S3 [H1],
confound, convergence), their figures, and a summary JSON, under
results/<dataset>/. Datasets without a reference annotation simply skip the
annotation-dependent analyses.

Usage (from project root):
    PYTHONPATH=. python scripts/analyze_dataset.py --dataset breast_cancer
"""

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.statistics import generate_stats_figures, run_all  # noqa: E402

logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="visium_hne",
                        help="dataset name; results read from results/ or results/<dataset>/")
    args = parser.parse_args()

    project = Path(__file__).resolve().parent.parent
    results_dir = project / "results" if args.dataset == "visium_hne" \
        else project / "results" / args.dataset

    if not (results_dir / "aggregate" / "per_spot_stability.csv").exists():
        logger.error("No sweep aggregates in %s; run the sweep first.", results_dir)
        sys.exit(1)

    summary = run_all(results_dir, results_dir / "aggregate" / "statistics_summary.json")
    generate_stats_figures(results_dir, results_dir / "figures")

    logger.info("=== %s summary ===", args.dataset)
    logger.info(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
