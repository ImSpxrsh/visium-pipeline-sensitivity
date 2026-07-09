"""
Generate all publication figures from a completed sweep.

Run from the project root:  PYTHONPATH=. python scripts/generate_figures.py
Requires that scripts/main.py (or scripts/aggregate.py) has populated results/.
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.visualization import generate_all_figures  # noqa: E402

logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project_dir = Path(__file__).resolve().parent.parent
    results_dir = project_dir / "results"
    generate_all_figures(results_dir, results_dir / "figures", verbose=True)
    logger.info("Figure generation complete")


if __name__ == "__main__":
    main()
