"""
Generate publication figures from a completed sweep.

Run from the project root:  python scripts/generate_figures.py
Requires that scripts/main.py has already populated results/.
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import load_json  # noqa: E402
from src.visualization import (  # noqa: E402
    generate_all_figures,
    plot_specification_curve,
)

logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project_dir = Path(__file__).resolve().parent.parent
    results_dir = project_dir / "results"
    figures_dir = results_dir / "figures"

    generate_all_figures(results_dir, figures_dir, verbose=True)

    # Specification curve from the persisted per-spec metadata.
    spec_curve_path = results_dir / "aggregate" / "specification_curves.json"
    if spec_curve_path.exists():
        metadata_list = load_json(spec_curve_path)["metadata"]
        plot_specification_curve(
            metadata_list,
            metric="n_clusters_found",
            output_path=figures_dir / "03_spec_curve_n_clusters.png",
        )
    else:
        logger.warning("No specification_curves.json found; run the sweep first")

    logger.info("Figure generation complete")


if __name__ == "__main__":
    main()
