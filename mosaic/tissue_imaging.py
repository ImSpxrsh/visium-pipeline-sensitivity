"""
Real histology imaging for MOSAIC figures.

Every Visium AnnData ships the actual H&E microscopy image the section was
photographed on. These helpers render that real tissue: the raw histology,
magnified boundary/core insets with spot overlays, and single-marker-gene
spatial expression (pseudo-microscopy). Coordinates in ``obsm['spatial']`` are
full-resolution pixels; the hires image is scaled by ``tissue_hires_scalef``.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

logger = logging.getLogger(__name__)

# Visium spot diameter is 55 um; spot_diameter_fullres gives it in full-res px.
SPOT_DIAMETER_UM = 55.0
S3_CMAP = "RdYlGn"


def _spatial_info(adata, library_id: Optional[str] = None) -> Dict[str, Any]:
    """Return the hires image, its scale factor, and um-per-fullres-pixel."""
    spatial = adata.uns["spatial"]
    library_id = library_id or next(iter(spatial))
    sp = spatial[library_id]
    hires = np.clip(np.asarray(sp["images"]["hires"]), 0, 1)
    scalef = sp["scalefactors"]["tissue_hires_scalef"]
    um_per_fullres_px = SPOT_DIAMETER_UM / sp["scalefactors"]["spot_diameter_fullres"]
    return {
        "image": hires,
        "scalef": scalef,
        "um_per_hires_px": um_per_fullres_px / scalef,
        "coords_hires": adata.obsm["spatial"] * scalef,
    }


def _add_scale_bar(ax, um_per_px: float, length_um: float = 500.0, loc=(0.05, 0.10)):
    """Draw a white scale bar of ``length_um`` at axes-fraction ``loc``."""
    length_px = length_um / um_per_px
    x0, y0 = ax.transLimits.inverted().transform(loc)
    ax.plot([x0, x0 + length_px], [y0, y0], color="white", lw=3, solid_capstyle="butt")
    ax.text(x0 + length_px / 2, y0, f"{length_um:.0f} um", color="white",
            ha="center", va="bottom", fontsize=8, fontweight="bold")


def plot_raw_histology(
    adata,
    output_path: Optional[Path] = None,
    library_id: Optional[str] = None,
    title: str = "H&E histology (raw tissue)",
) -> None:
    """Show the full-resolution H&E image with no analytical overlay."""
    info = _spatial_info(adata, library_id)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(info["image"])
    _add_scale_bar(ax, info["um_per_hires_px"])
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axis("off")
    _save(fig, output_path)


def _pick_core_and_boundary(per_spot: pd.DataFrame, coords_hires: np.ndarray,
                            obs_names, crop_px: int = 220) -> Dict[str, Tuple[float, float]]:
    """
    Choose a robust (high-S3) and a fragile (low-S3) interior location.

    Uses a locally smoothed S3 (mean over ~30 nearest spots) so the picked
    centers sit in coherent zones, and requires the two centers to be at least
    two crop widths apart so their call-outs do not overlap.
    """
    from scipy.spatial import cKDTree

    s3 = per_spot.set_index("spot_id").reindex(obs_names)["s3_score"].values
    tree = cKDTree(coords_hires)
    _, nbr = tree.query(coords_hires, k=min(30, len(coords_hires)))
    local = np.nanmean(s3[nbr], axis=1)

    x, y = coords_hires[:, 0], coords_hires[:, 1]
    interior = (
        (x > np.percentile(x, 12)) & (x < np.percentile(x, 88))
        & (y > np.percentile(y, 12)) & (y < np.percentile(y, 88))
    )
    idx = np.where(interior)[0]
    robust = idx[np.argmax(local[idx])]
    rc = coords_hires[robust]
    # Fragile center: lowest local S3 among interior spots far from the robust one.
    far = idx[np.linalg.norm(coords_hires[idx] - rc, axis=1) > 2.2 * crop_px]
    pool = far if len(far) else idx
    fragile = pool[np.argmin(local[pool])]
    return {
        "robust core": (coords_hires[robust, 0], coords_hires[robust, 1]),
        "fragile boundary": (coords_hires[fragile, 0], coords_hires[fragile, 1]),
    }


def plot_zoomed_inset(
    adata,
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
    library_id: Optional[str] = None,
    crop_px: int = 220,
) -> None:
    """
    Magnified H&E crops at a robust core and a fragile boundary, with S3 overlays.

    Left: full tissue with dashed call-out boxes. Right: the two crops, each with
    spots colored by S3, showing that fragile zones are morphological transitions.
    """
    info = _spatial_info(adata, library_id)
    img, coords = info["image"], info["coords_hires"]
    s3 = per_spot.set_index("spot_id").reindex(adata.obs_names)["s3_score"].values
    regions = _pick_core_and_boundary(per_spot, coords, adata.obs_names, crop_px)

    fig = plt.figure(figsize=(15, 7))
    gsL = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.15)
    axFull = fig.add_subplot(gsL[0])
    axFull.imshow(img)
    axFull.axis("off")
    axFull.set_title("A  Full section with call-outs", fontsize=11, fontweight="bold")
    _add_scale_bar(axFull, info["um_per_hires_px"])

    gsR = gsL[1].subgridspec(2, 1, hspace=0.2)
    for row, (label, (cx, cy)) in enumerate(regions.items()):
        axFull.add_patch(Rectangle((cx - crop_px, cy - crop_px), 2 * crop_px, 2 * crop_px,
                                   fill=False, edgecolor="white", lw=2, linestyle="--"))
        axFull.text(cx, cy - crop_px - 10, label, color="white", fontsize=9,
                    ha="center", va="bottom", fontweight="bold")

        ax = fig.add_subplot(gsR[row])
        ax.imshow(img)
        m = ((coords[:, 0] > cx - crop_px) & (coords[:, 0] < cx + crop_px)
             & (coords[:, 1] > cy - crop_px) & (coords[:, 1] < cy + crop_px))
        sc = ax.scatter(coords[m, 0], coords[m, 1], c=s3[m], cmap=S3_CMAP,
                        vmin=0, vmax=1, s=90, edgecolors="white", linewidths=0.5)
        ax.set_xlim(cx - crop_px, cx + crop_px)
        ax.set_ylim(cy + crop_px, cy - crop_px)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{'BC'[row]}  {label}", fontsize=11, fontweight="bold")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label("S3", fontsize=8)

    fig.suptitle("Real tissue at fragile boundaries vs robust cores",
                 fontsize=13, fontweight="bold")
    _save(fig, output_path)


def plot_marker_gene_expression(
    adata,
    genes: List[str],
    output_path: Optional[Path] = None,
    library_id: Optional[str] = None,
) -> None:
    """Spatial expression of marker genes on the H&E image (pseudo-microscopy)."""
    import scanpy as sc

    info = _spatial_info(adata, library_id)
    img, coords = info["image"], info["coords_hires"]

    a = adata.copy()
    if "counts" in a.layers:
        a.X = a.layers["counts"].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    genes = [g for g in genes if g in a.var_names]
    n = len(genes)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for i, g in enumerate(genes):
        ax = axes.ravel()[i]
        ax.axis("off")
        ax.imshow(img)
        col = a[:, g].X
        expr = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
        sc_ = ax.scatter(coords[:, 0], coords[:, 1], c=expr, cmap="magma", s=10)
        ax.set_title(g, fontsize=11, fontweight="bold", style="italic")
        fig.colorbar(sc_, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Marker gene spatial expression", fontsize=13, fontweight="bold")
    _save(fig, output_path)


def _save(fig, output_path: Optional[Path]) -> None:
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)
