"""
Load and validate spatial transcriptomics datasets.

Primary dataset: squidpy's Visium H&E demo, a 10x Genomics **mouse brain**
coronal section (~2,688 spots x ~18,078 genes) with an associated H&E image.
It ships with a reference anatomical annotation in ``adata.obs['cluster']``
(Cortex, Hippocampus, Thalamus, Striatum, Fiber tract, ...), which we use as
ground truth for validation.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)


def load_dataset(
    dataset_name: str = "visium_hne",
    verbose: bool = True,
) -> Tuple[ad.AnnData, Optional[np.ndarray]]:
    """
    Load a public spatial transcriptomics dataset via squidpy.

    Args:
        dataset_name: "visium_hne" (mouse brain, H&E) or
                      "visium_fluo" (mouse brain, fluorescence).
        verbose: Print progress.

    Returns:
        (adata, image) where ``image`` is the H&E array if available, else None.

    Raises:
        ValueError: If the dataset name is not recognized.
    """
    import squidpy as sq  # imported lazily so unit tests can run without it

    if dataset_name == "visium_hne":
        if verbose:
            logger.info("Loading Visium H&E (mouse brain, coronal)...")
        adata = sq.datasets.visium_hne_adata()
        # The H&E image (~400 MB) is used only for figure overlays, not the
        # analysis. Fetch it best-effort so a slow/failed download never blocks
        # the multiverse sweep.
        image = _try_load_image(sq.datasets.visium_hne_image, verbose=verbose)
    elif dataset_name == "visium_fluo":
        if verbose:
            logger.info("Loading Visium fluorescence (mouse brain)...")
        adata = sq.datasets.visium_fluo_adata()
        image = None
    elif dataset_name in TENX_VISIUM or dataset_name in TENX_ALIASES:
        sample_id = TENX_ALIASES.get(dataset_name, dataset_name)
        if verbose:
            logger.info("Loading 10x Visium %s ...", sample_id)
        # Public 10x Visium (human tumor/disease). Raw filtered_feature_bc_matrix;
        # the pipeline re-derives normalization from counts, preserving Year-1
        # invariants. The tissue image is embedded in adata.uns['spatial'].
        adata = sq.datasets.visium(sample_id)
        adata.var_names_make_unique()
        image = None
    else:
        raise ValueError(
            f"Unknown dataset: {dataset_name!r}. Known: visium_hne, visium_fluo, "
            f"aliases {sorted(TENX_ALIASES)}, or any 10x id in {sorted(TENX_VISIUM)}."
        )

    return adata, image


# Public 10x Visium human tumor/disease datasets reachable without registration.
TENX_VISIUM = {
    "V1_Breast_Cancer_Block_A_Section_1",
    "V1_Breast_Cancer_Block_A_Section_2",
    "Parent_Visium_Human_Glioblastoma",
    "Targeted_Visium_Human_Glioblastoma_Pan_Cancer",
}
# Friendly aliases so results namespace cleanly (results/<alias>/...).
TENX_ALIASES = {
    "breast_cancer": "V1_Breast_Cancer_Block_A_Section_1",
    "breast_cancer_2": "V1_Breast_Cancer_Block_A_Section_2",
    "glioblastoma": "Parent_Visium_Human_Glioblastoma",
}


def _try_load_image(loader, retries: int = 2, verbose: bool = True):
    """Attempt to load the tissue image, tolerating download failures."""
    for attempt in range(1, retries + 1):
        try:
            return loader()
        except Exception as e:  # network timeout, partial download, etc.
            if verbose:
                logger.warning(
                    "Image load failed (attempt %d/%d): %s", attempt, retries, e
                )
    logger.warning("Proceeding without the H&E image; overlays will use spot coords.")
    return None


def validate_dataset(adata: ad.AnnData, verbose: bool = True) -> bool:
    """
    Check that an AnnData object has the fields the pipeline requires.

    Returns True only if every check passes.
    """
    checks = {
        "Has expression matrix (X)": adata.X is not None,
        "Has spatial coordinates": "spatial" in adata.obsm,
        "Has feature names": adata.var_names is not None and adata.n_vars > 0,
        "Has observation names": adata.obs_names is not None and adata.n_obs > 0,
        "Positive spot count": adata.n_obs > 0,
        "Positive gene count": adata.n_vars > 0,
    }

    all_pass = all(checks.values())

    if verbose:
        logger.info("Dataset validation (%d spots x %d genes):", adata.n_obs, adata.n_vars)
        for check, result in checks.items():
            logger.info("  %s %s", "PASS" if result else "FAIL", check)

    return all_pass


def prepare_data(
    adata: ad.AnnData,
    output_dir: Optional[Path] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Minimal, leakage-free preparation shared by every spec.

    Ensures unique gene names and stashes the raw counts in
    ``adata.layers["counts"]`` so each spec can normalize independently.
    Does NOT normalize or filter on global statistics.

    Args:
        adata: Raw AnnData object.
        output_dir: If given, write the prepared object as an .h5ad here.
        verbose: Print progress.

    Returns:
        A prepared copy of ``adata``.
    """
    adata = adata.copy()

    if verbose:
        logger.info("Preparing data: %d spots x %d genes", adata.n_obs, adata.n_vars)

    adata.var_names_make_unique()

    # Preserve raw counts; each spec re-derives X from this layer.
    adata.layers["counts"] = adata.X.copy() if hasattr(adata.X, "copy") else adata.X

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "visium_prepared.h5ad"
        adata.write_h5ad(output_path)
        if verbose:
            logger.info("Saved prepared data to %s", output_path)

    return adata


def main() -> None:
    """Smoke-test: load, validate, and prepare the Visium H&E demo."""
    logging.basicConfig(level=logging.INFO)

    adata, image = load_dataset("visium_hne", verbose=True)
    logger.info("Loaded: %d spots x %d genes", adata.n_obs, adata.n_vars)
    logger.info("H&E image present: %s", image is not None)

    if not validate_dataset(adata, verbose=True):
        raise ValueError("Dataset validation failed")

    prepare_data(adata, output_dir=Path("data/processed"), verbose=True)
    logger.info("Data loading complete")


if __name__ == "__main__":
    main()
