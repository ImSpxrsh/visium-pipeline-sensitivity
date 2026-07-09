"""
Generate the real-histology figures (raw H&E, boundary insets, marker genes).

Run from the project root:  PYTHONPATH=. python scripts/generate_tissue.py
Requires the primary sweep aggregates (results/aggregate/per_spot_stability.csv).
"""

import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mosaic.data_loader import load_dataset, prepare_data  # noqa: E402
from mosaic.tissue_imaging import (  # noqa: E402
    plot_marker_gene_expression,
    plot_raw_histology,
    plot_zoomed_inset,
)

logger = logging.getLogger(__name__)

# Canonical mouse-brain region markers (fiber tracts, hippocampus, thalamus, cortex).
MARKERS = ["Snap25", "Mbp", "Plp1", "Hpca", "Pcp4", "Nrgn"]


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    figures = project / "results" / "figures"

    adata, _ = load_dataset("visium_hne", verbose=False)
    adata = prepare_data(adata, verbose=False)
    per_spot = pd.read_csv(project / "results" / "aggregate" / "per_spot_stability.csv")

    plot_raw_histology(adata, figures / "01_raw_histology.png")
    plot_zoomed_inset(adata, per_spot, figures / "05_zoomed_insets.png")
    present = [g for g in MARKERS if g in adata.var_names]
    plot_marker_gene_expression(adata, present, figures / "06_marker_genes.png")
    logger.info("Tissue figures written.")


if __name__ == "__main__":
    main()
