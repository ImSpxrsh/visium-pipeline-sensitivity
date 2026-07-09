"""
Generate the cross-dataset generalization figure (Extended Data 8).

Requires that scripts/main.py has been run for both the primary dataset
(results/) and at least one secondary dataset (results/<dataset>/).

Run from the project root:
    PYTHONPATH=. python scripts/generate_cross_dataset.py
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.visualization import ed8_cross_dataset  # noqa: E402

logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    results = project / "results"

    datasets = [{"label": "Mouse brain (H&E)", "results_dir": results}]
    fluo = results / "visium_fluo"
    if (fluo / "aggregate" / "per_spot_stability.csv").exists():
        datasets.append({"label": "Mouse brain (fluorescence)", "results_dir": fluo})
    else:
        logger.warning("No visium_fluo results found; ED8 will show one dataset only.")

    ed8_cross_dataset(datasets, output_path=results / "figures" / "ed8_cross_dataset.png")
    logger.info("Cross-dataset figure written.")


if __name__ == "__main__":
    main()
