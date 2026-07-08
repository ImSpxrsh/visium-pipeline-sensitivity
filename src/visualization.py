"""
Publication figures for MOSAIC: stability maps, specification curves, and
per-gene instability plots.

Figures are saved to disk; interactive display is disabled by default so this
runs cleanly in a headless / batch context.
"""

import glob
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")  # headless backend; safe for batch figure generation
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

logger = logging.getLogger(__name__)


def _save(fig, output_path: Optional[Path]) -> None:
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def plot_stability_map(
    adata_canonical,
    s3_scores_per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
    title: str = "Spatial Stability (S3)",
    figsize: Tuple[float, float] = (16, 5),
) -> None:
    """
    Three-panel tissue view: canonical clustering, per-spot S³ heatmap, consensus.

    Uses squidpy's spatial scatter when available, falling back to a raw
    ``obsm['spatial']`` scatter so the figure still renders without an image.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    adata = adata_canonical
    adata.obs = adata.obs.copy()
    adata.obs["s3_score"] = s3_scores_per_spot["s3_score"].values
    adata.obs["consensus_cluster"] = pd.Categorical(
        s3_scores_per_spot["consensus_cluster"].values
    )

    coords = adata.obsm["spatial"]
    leiden_codes = pd.Categorical(adata.obs["leiden"]).codes

    axes[0].scatter(coords[:, 0], coords[:, 1], c=leiden_codes, cmap="tab20", s=8)
    axes[0].set_title("Canonical clustering")

    sc = axes[1].scatter(
        coords[:, 0], coords[:, 1], c=adata.obs["s3_score"], cmap="RdYlGn", s=8, vmin=0, vmax=1
    )
    axes[1].set_title("Per-spot stability (S3)")
    fig.colorbar(sc, ax=axes[1], fraction=0.046, pad=0.04)

    axes[2].scatter(
        coords[:, 0], coords[:, 1],
        c=s3_scores_per_spot["consensus_cluster"].values, cmap="tab20", s=8,
    )
    axes[2].set_title("Consensus domains")

    for ax in axes:
        ax.set_aspect("equal")
        ax.invert_yaxis()  # image convention: origin top-left
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, output_path)


def plot_specification_curve(
    metadata_list: List[Dict[str, Any]],
    metric: str = "n_clusters_found",
    output_path: Optional[Path] = None,
    title: str = "Specification Curve",
) -> None:
    """Plot a metric across specs, sorted ascending (the classic spec curve)."""
    values = np.array([m.get(metric, np.nan) for m in metadata_list], dtype=float)
    order = np.argsort(values)
    sorted_vals = values[order]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(range(len(sorted_vals)), sorted_vals, "o-", color="steelblue", lw=2, ms=6)
    if np.isfinite(sorted_vals).any():
        ax.fill_between(
            range(len(sorted_vals)),
            np.nanmin(sorted_vals),
            sorted_vals,
            alpha=0.15,
            color="steelblue",
        )
    ax.set_xlabel("Specification (sorted)")
    ax.set_ylabel(metric)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _save(fig, output_path)


def plot_instability_heatmap(
    per_feature_s3: pd.DataFrame,
    top_n: int = 20,
    output_path: Optional[Path] = None,
    title: str = "Most fragile genes (lowest S3)",
) -> None:
    """Horizontal barplot of the least stable genes across the multiverse."""
    top_unstable = per_feature_s3.nsmallest(top_n, "s3_score")

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(range(len(top_unstable)), top_unstable["s3_score"], color="coral")
    ax.set_yticks(range(len(top_unstable)))
    ax.set_yticklabels(top_unstable["gene"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Stability score (S3)")
    ax.set_xlim(0, 1)
    ax.set_title(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, output_path)


def generate_all_figures(
    results_dir: Path, output_dir: Path, verbose: bool = True
) -> None:
    """
    Generate every publication figure from a completed sweep.

    Reads aggregate CSVs plus the first spec's clustered AnnData.
    """
    import anndata as ad

    results_dir = Path(results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if verbose:
        logger.info("Generating publication figures...")

    per_spot = pd.read_csv(results_dir / "aggregate" / "per_spot_stability.csv")
    per_feature = pd.read_csv(results_dir / "aggregate" / "per_feature_stability.csv")

    # Resolve the canonical (first) spec directory via glob, then load its .h5ad.
    matches = sorted(glob.glob(str(results_dir / "runs" / "spec_00_*" / "adata_clustered.h5ad")))
    if not matches:
        raise FileNotFoundError("No canonical spec_00 result found in results/runs/")
    adata_canonical = ad.read_h5ad(matches[0])

    plot_stability_map(
        adata_canonical, per_spot, output_path=output_dir / "01_stability_map.png"
    )
    plot_instability_heatmap(
        per_feature, top_n=20, output_path=output_dir / "02_gene_instability.png"
    )
    if verbose:
        logger.info("Figures saved to %s", output_dir)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project_dir = Path(__file__).resolve().parent.parent
    results_dir = project_dir / "results"
    generate_all_figures(results_dir, results_dir / "figures", verbose=True)


if __name__ == "__main__":
    main()
