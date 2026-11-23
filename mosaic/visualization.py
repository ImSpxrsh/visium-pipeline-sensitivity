"""
Publication figures for MOSAIC.

Produces the paper's main and extended-data figures from a completed sweep
(results/runs/spec_*/ plus results/aggregate/*). All static figures use a
headless Agg backend and colorblind-safe palettes; the interactive multiverse
explorer is written as a self-contained Plotly HTML file.

Figures that would require data outside this study (pathologist ground-truth
annotations, additional tissues/modalities) are intentionally not fabricated.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

logger = logging.getLogger(__name__)

# --- palettes (colorblind-safe) ------------------------------------------------
# Paul Tol / ColorBrewer qualitative for categorical cluster labels.
CLUSTER_COLORS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFD92F",
    "#A65628", "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3",
    "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7",
    "#BEBADA", "#FB8072",
]
CLUSTER_CMAP = ListedColormap(CLUSTER_COLORS)
S3_CMAP = "RdYlGn"  # red (fragile) -> yellow -> green (robust)

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


# --- data loading --------------------------------------------------------------
def _slim_adata(adata):
    import anndata as ad

    # Drop uns entirely: it holds the tissue images (~50 MB per spec), and the
    # analyses read SVG tables from the per-spec CSV (run["svg_df"]) instead.
    return ad.AnnData(obs=adata.obs.copy(), obsm=dict(adata.obsm.items()))


def load_runs(results_dir: Path, lean: bool = False) -> List[Dict[str, Any]]:
    """
    Load every completed spec directory into a list of run records.

    Each record: {idx, name, spec, adata, svg_df, metadata}. Records are sorted
    by spec index. The (small) SVG table is read from the per-spec CSV.

    With ``lean=True``, every run except the first keeps only obs, obsm, and uns
    (the analyses read full data only from ``runs[0]``), so sweeps with hundreds
    of specs fit in memory.
    """
    import anndata as ad

    results_dir = Path(results_dir)
    runs = []
    for spec_dir in sorted((results_dir / "runs").glob("spec_*")):
        adata_path = spec_dir / "adata_clustered.h5ad"
        meta_path = spec_dir / "metadata.json"
        if not (adata_path.exists() and meta_path.exists()):
            continue
        metadata = json.loads(meta_path.read_text())
        svg_csv = spec_dir / "svg_results.csv"
        svg_df = pd.read_csv(svg_csv, index_col=0) if svg_csv.exists() else None
        idx = int(spec_dir.name.split("_")[1])
        runs.append(
            {
                "idx": idx,
                "name": metadata.get("spec_name", spec_dir.name),
                "spec": metadata.get("spec", {}),
                "adata": _slim_adata(ad.read_h5ad(adata_path)) if lean and runs
                else ad.read_h5ad(adata_path),
                "svg_df": svg_df,
                "metadata": metadata,
            }
        )
    logger.info("Loaded %d completed runs", len(runs))
    return runs


def _by_name(runs: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    return {r["name"]: r for r in runs}


def _significant(svg_df: pd.DataFrame, alpha: float = 0.05) -> pd.Series:
    """Boolean SVG calls indexed by gene, using the best available p-value column."""
    for col in ("pval_norm_fdr_bh", "pval_sim_fdr_bh", "pval_norm", "pval_sim", "pval"):
        if col in svg_df.columns:
            return pd.Series((svg_df[col] < alpha).values, index=svg_df.index)
    return pd.Series((svg_df["I"] > svg_df["I"].median()).values, index=svg_df.index)


def svg_concordance_matrix(
    runs: List[Dict[str, Any]], alpha: float = 0.05, common_only: bool = False
) -> pd.DataFrame:
    """
    Gene x spec boolean matrix: is each gene a significant SVG in each spec?

    Specs use different HVG sets, so a gene may be untested in some specs. With
    ``common_only=True`` we keep only genes tested in *every* spec, isolating
    genuine SVG-test flips from feature-selection differences (a gene absent from
    a spec's HVG set is otherwise drawn as "not called"). With ``common_only=False``
    absent genes are treated as not-called, matching the full-multiverse S3.
    """
    cols = {}
    for r in runs:
        if r["svg_df"] is None:
            continue
        cols[r["name"]] = _significant(r["svg_df"], alpha)
    mat = pd.DataFrame(cols)  # union of genes; NaN where a gene was not tested
    if common_only:
        mat = mat.dropna(axis=0, how="any")
    return mat.fillna(False).astype(bool)


# --- shared spatial-scatter helper --------------------------------------------
def _spatial_scatter(ax, coords, values, *, categorical, title, s=6, cmap=None, vmin=None, vmax=None):
    if categorical:
        codes = pd.Categorical(values).codes
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=codes % len(CLUSTER_COLORS),
                        cmap=CLUSTER_CMAP, s=s, vmin=0, vmax=len(CLUSTER_COLORS) - 1)
    else:
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, cmap=cmap or S3_CMAP,
                        s=s, vmin=vmin, vmax=vmax)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title)
    return sc


# ==============================================================================
# Figure 1 -- 6-panel tissue montage
# ==============================================================================
def figure1_spec_montage(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> None:
    """Five specification clusterings + the per-spot S3 overlay (2x3 grid)."""
    named = _by_name(runs)
    wanted = [
        ("canonical_seurat", "Canonical: Seurat 3k HVGs, Leiden res=0.5"),
        ("conservative_2k_hvgs", "Conservative: 2k HVGs"),
        ("aggressive_5k_hvgs", "Aggressive: 5k HVGs"),
        ("coarse_res03", "Coarse: Leiden res=0.3"),
        ("fine_res08", "Fine: Leiden res=0.8"),
    ]
    # Fall back to the first runs if a named spec is missing.
    picks = [named.get(n) for n, _ in wanted]
    picks = [p if p is not None else runs[i] for i, p in enumerate(picks)]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    for ax, (run, (_, label)) in zip(axes[:5], zip(picks, wanted)):
        a = run["adata"]
        n = int(a.obs["leiden"].nunique())
        _spatial_scatter(ax, a.obsm["spatial"], a.obs["leiden"].values,
                         categorical=True, title=f"{label}\n-> {n} domains")
        ax.title.set_fontsize(10)

    # Panel F: S3 overlay on the canonical geometry.
    ref = picks[0]["adata"]
    s3 = per_spot.set_index("spot_id").reindex(ref.obs_names)["s3_score"].values
    scF = _spatial_scatter(axes[5], ref.obsm["spatial"], s3, categorical=False,
                           title="Spatial Stability Score (S3)", cmap=S3_CMAP, vmin=0, vmax=1)
    cb = fig.colorbar(scF, ax=axes[5], fraction=0.046, pad=0.04)
    cb.set_label("S3 (green robust, red fragile)")
    axes[5].title.set_fontsize(10)

    for i, ax in enumerate(axes):
        ax.text(-0.08, 1.14, chr(ord("A") + i), transform=ax.transAxes,
                fontsize=16, fontweight="bold", va="top", ha="left")

    fig.suptitle(
        "Biological conclusions depend sensitively on analytical choices",
        fontsize=15, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, output_path)


# ==============================================================================
# Figure 2 -- specification-curve matrix
# ==============================================================================
def figure2_spec_curve_matrix(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
    top_genes: int = 25,
) -> None:
    """Cluster count, spec-consensus agreement, and SVG concordance across specs."""
    from mosaic.stability import _align_labels_to_reference

    order = sorted(range(len(runs)), key=lambda i: runs[i]["metadata"].get("n_clusters_found", 0))
    runs_ord = [runs[i] for i in order]
    names = [r["name"] for r in runs_ord]
    n_clusters = [r["metadata"].get("n_clusters_found", np.nan) for r in runs_ord]

    # Row 2: each spec's agreement with the multiverse consensus (label-aligned).
    ref = runs_ord[0]["adata"]
    consensus = per_spot.set_index("spot_id").reindex(ref.obs_names)["consensus_cluster"].values
    agreement = []
    for r in runs_ord:
        a = r["adata"]
        if not a.obs_names.equals(ref.obs_names):
            agreement.append(np.nan)
            continue
        aligned = _align_labels_to_reference(a.obs["leiden"].values, ref.obs["leiden"].values)
        agreement.append(float(np.mean(aligned == consensus)))

    # Row 3: SVG concordance for genes tested in every spec (isolates SVG-test
    # flips from feature-selection differences between specs).
    conc = svg_concordance_matrix(runs_ord, common_only=True)
    if not conc.empty:
        churn = (conc.mean(axis=1) - 0.5).abs()  # genes that flip most sit near 0.5
        genes = churn.nsmallest(top_genes).index
        conc_top = conc.loc[genes, names] if set(names).issubset(conc.columns) else conc.loc[genes]
    else:
        conc_top = pd.DataFrame()

    fig = plt.figure(figsize=(max(10, 0.7 * len(runs_ord)), 11))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.1, 0.35, 2.2], hspace=0.35)
    x = np.arange(len(runs_ord))

    ax1 = fig.add_subplot(gs[0])
    ax1.bar(x, n_clusters, color="#377EB8")
    ax1.set_ylabel("N domains")
    ax1.set_title("Headline conclusions across the analytical multiverse", loc="left")
    ax1.set_xticks([])

    ax2 = fig.add_subplot(gs[1])
    ax2.imshow(np.array(agreement)[None, :], aspect="auto", cmap=S3_CMAP, vmin=0, vmax=1,
               extent=[-0.5, len(runs_ord) - 0.5, 0, 1])
    ax2.set_yticks([0.5])
    ax2.set_yticklabels(["consensus\nagreement"])
    ax2.set_xticks([])

    ax3 = fig.add_subplot(gs[2])
    if not conc_top.empty:
        ax3.imshow(conc_top.values, aspect="auto", cmap="binary", interpolation="nearest")
        ax3.set_yticks(range(len(conc_top)))
        ax3.set_yticklabels(conc_top.index, fontsize=7)
        ax3.set_ylabel("Spatially variable gene")
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{r['idx']}:{r['name'][:12]}" for r in runs_ord],
                        rotation=90, fontsize=7)
    ax3.set_title("SVG called (black) vs not called (white) per specification", loc="left")

    fig.suptitle("Figure 2 | Specification-curve matrix", fontsize=14, fontweight="bold")
    _save(fig, output_path)


# ==============================================================================
# Figure 3 -- stability map + summary statistics
# ==============================================================================
def figure3_stability_summary(
    canonical_adata,
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> None:
    """Large S3 map, S3 histogram, and mean S3 by data-driven consensus domain."""
    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.4, 1, 1], wspace=0.35)

    s3 = per_spot.set_index("spot_id").reindex(canonical_adata.obs_names)
    coords = canonical_adata.obsm["spatial"]

    axA = fig.add_subplot(gs[0])
    scA = _spatial_scatter(axA, coords, s3["s3_score"].values, categorical=False,
                           title="A  Per-spot S3", cmap=S3_CMAP, vmin=0, vmax=1, s=10)
    fig.colorbar(scA, ax=axA, fraction=0.046, pad=0.04).set_label("S3")

    axB = fig.add_subplot(gs[1])
    vals = per_spot["s3_score"].values
    axB.hist(vals, bins=25, color="#999999", edgecolor="white")
    for thr, c in [(0.5, "#D73027"), (0.8, "#1A9850")]:
        axB.axvline(thr, color=c, ls="--", lw=1.5)
    robust = float(np.mean(vals >= 0.8)) * 100
    fragile = float(np.mean(vals < 0.5)) * 100
    axB.set_title(f"B  S3 distribution\n{robust:.0f}% robust, {fragile:.0f}% fragile")
    axB.set_xlabel("S3")
    axB.set_ylabel("Spots")

    axC = fig.add_subplot(gs[2])
    dom = per_spot.groupby("consensus_cluster")["s3_score"].agg(["mean", "std", "count"])
    dom = dom.sort_values("mean")
    axC.barh(range(len(dom)), dom["mean"], xerr=dom["std"].fillna(0),
             color="#4DAF4A", ecolor="#555555", capsize=2)
    axC.set_yticks(range(len(dom)))
    axC.set_yticklabels([f"domain {int(i)} (n={int(c)})" for i, c in zip(dom.index, dom["count"])],
                        fontsize=8)
    axC.set_xlim(0, 1)
    axC.set_xlabel("Mean S3")
    axC.set_title("C  Stability by consensus domain")

    fig.suptitle("Figure 3 | Spatial Stability Score reveals robust vs contingent structure",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ==============================================================================
# Extended Data 4 -- specification space table
# ==============================================================================
def ed4_spec_space_table(
    runs: List[Dict[str, Any]],
    output_path: Optional[Path] = None,
) -> None:
    """Render the curated specification space as a table figure."""
    rows = []
    for r in runs:
        s = r["spec"]
        rows.append([
            r["idx"], r["name"],
            s.get("normalization", ""),
            f"{s.get('feature_selection', {}).get('method', '')}/"
            f"{s.get('feature_selection', {}).get('n_genes', '')}",
            s.get("dimensionality_reduction", {}).get("n_components", ""),
            f"{s.get('neighborhood_graph', {}).get('type', '')} "
            f"k={s.get('neighborhood_graph', {}).get('k_or_radius', '')}",
            f"{s.get('clustering', {}).get('method', '')} "
            f"res={s.get('clustering', {}).get('resolution', '')}",
            r["metadata"].get("n_clusters_found", ""),
        ])
    cols = ["#", "name", "norm", "HVG", "PCA", "graph", "clustering", "domains"]

    fig, ax = plt.subplots(figsize=(13, 0.32 * len(rows) + 1.2))
    ax.axis("off")
    # Give the name column extra room; keep the rest compact.
    col_widths = [0.05, 0.20, 0.09, 0.13, 0.07, 0.13, 0.16, 0.09]
    tbl = ax.table(cellText=rows, colLabels=cols, loc="upper center", cellLoc="left",
                   colWidths=col_widths)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.25)
    for j in range(len(cols)):
        tbl[0, j].set_facecolor("#377EB8")
        tbl[0, j].set_text_props(color="white", fontweight="bold")
    ax.set_title("Extended Data 4 | Curated specification space",
                 fontsize=13, fontweight="bold", pad=12)
    _save(fig, output_path)


# ==============================================================================
# Extended Data 5 -- per-feature (SVG) stability
# ==============================================================================
def ed5_per_feature_stability(
    per_feature: pd.DataFrame,
    runs: List[Dict[str, Any]],
    output_path: Optional[Path] = None,
    top_n: int = 20,
) -> None:
    """Most-contested SVGs barplot + gene x spec SVG concordance heatmap.

    Restricted to genes tested in every specification, so "instability" reflects
    genuine significance flips rather than a gene being dropped by feature
    selection. Genes are ranked by how close their call frequency is to 50%
    (maximally contested); the SVG stability shown is that call frequency.
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={"width_ratios": [1, 1.3]})

    conc = svg_concordance_matrix(runs, common_only=True)
    if conc.empty:
        logger.warning("No genes tested in every spec; skipping ED5")
        plt.close(fig)
        return

    gene_freq = conc.mean(axis=1)  # fraction of specs calling each gene an SVG
    contested = (gene_freq - 0.5).abs().nsmallest(top_n).index
    freq = gene_freq.loc[contested].sort_values()

    colors = ["#D73027" if v < 0.5 else "#FEE08B" if v < 0.8 else "#1A9850" for v in freq]
    axes[0].barh(range(len(freq)), freq.values, color=colors)
    axes[0].set_yticks(range(len(freq)))
    axes[0].set_yticklabels(freq.index, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1)
    axes[0].axvline(0.5, color="#555555", ls="--", lw=1)
    axes[0].set_xlabel("Fraction of specs calling gene an SVG")
    axes[0].set_title("A  Most contested SVG calls")

    sub = conc.loc[freq.index]
    axes[1].imshow(sub.values, aspect="auto", cmap="binary", interpolation="nearest")
    axes[1].set_yticks(range(len(sub)))
    axes[1].set_yticklabels(sub.index, fontsize=8)
    axes[1].set_xticks(range(sub.shape[1]))
    axes[1].set_xticklabels(sub.columns, rotation=90, fontsize=6)
    axes[1].set_title("B  SVG called (black) per specification")

    fig.suptitle("Extended Data 5 | Spatially variable gene calls are pipeline-contingent",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ==============================================================================
# Extended Data 6 -- headline instability metrics
# ==============================================================================
def ed6_headline_instability(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> None:
    """Cluster-count spread, top-SVG list overlap, and runtime vs stability."""
    from itertools import combinations

    n_clusters = np.array([r["metadata"].get("n_clusters_found", np.nan) for r in runs], float)
    resolutions = np.array(
        [r["spec"].get("clustering", {}).get("resolution", np.nan) for r in runs], float
    )
    methods = [r["spec"].get("clustering", {}).get("method", "?") for r in runs]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].boxplot(n_clusters[~np.isnan(n_clusters)], vert=True, widths=0.5)
    axes[0].scatter(np.random.RandomState(0).normal(1, 0.04, len(n_clusters)), n_clusters,
                    color="#377EB8", alpha=0.7, zorder=3)
    axes[0].set_xticks([])
    axes[0].set_ylabel("N domains")
    axes[0].set_title(f"A  Domain count\nmean {np.nanmean(n_clusters):.1f}, "
                      f"range {np.nanmin(n_clusters):.0f}-{np.nanmax(n_clusters):.0f}")

    # Top-10 SVG list overlap (Jaccard) between spec pairs.
    top_lists = []
    for r in runs:
        if r["svg_df"] is not None:
            top_lists.append(set(r["svg_df"].sort_values("I", ascending=False).head(10).index))
    jac = [len(a & b) / len(a | b) for a, b in combinations(top_lists, 2)] if len(top_lists) > 1 else [0]
    axes[1].hist(jac, bins=15, color="#984EA3", edgecolor="white")
    axes[1].axvline(np.mean(jac), color="k", ls="--")
    axes[1].set_xlabel("Jaccard overlap of top-10 SVGs (spec pairs)")
    axes[1].set_ylabel("Spec pairs")
    axes[1].set_title(f"B  Top-SVG list agreement\nmean Jaccard {np.mean(jac):.2f}")

    # Clustering resolution is the dominant lever on domain count; show it split by
    # algorithm. (Runtime is not shown: overnight machine-sleep makes the recorded
    # wall-clock unreliable as a compute measure.)
    method_colors = {"leiden": "#377EB8", "louvain": "#E41A1C"}
    for meth in sorted(set(methods)):
        idx = [i for i, m in enumerate(methods) if m == meth]
        jitter = np.random.RandomState(1).normal(0, 0.01, len(idx))
        axes[2].scatter(resolutions[idx] + jitter, n_clusters[idx],
                        c=method_colors.get(meth, "#666666"), s=60, alpha=0.8, label=meth)
    axes[2].set_xlabel("Clustering resolution")
    axes[2].set_ylabel("N domains")
    axes[2].set_title("C  Domain count vs resolution\n(color = algorithm)")
    axes[2].legend(fontsize=9)

    fig.suptitle("Extended Data 6 | Instability in headline conclusions",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, output_path)


# ==============================================================================
# Extended Data 7 -- validation against the reference anatomical annotation
# ==============================================================================
def ed7_reference_validation(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    output_path: Optional[Path] = None,
    reference_key: str = "cluster",
) -> None:
    """
    Validate multiverse output against the dataset's anatomical annotation.

    Panel A: normalized agreement between reference regions (rows) and consensus
    domains (cols). Panel B: mean per-spot S3 within each anatomical region,
    showing which known structures are recovered robustly.
    """
    ref_adata = runs[0]["adata"]
    if reference_key not in ref_adata.obs:
        logger.warning("No reference annotation %r; skipping ED7", reference_key)
        return

    ref = ref_adata.obs[reference_key].astype(str).values
    s3 = per_spot.set_index("spot_id").reindex(ref_adata.obs_names)
    consensus = s3["consensus_cluster"].astype(int).values

    # Contingency of reference region vs consensus domain, row-normalized.
    ct = pd.crosstab(pd.Series(ref, name="anatomical region"),
                     pd.Series(consensus, name="consensus domain"))
    ct_norm = ct.div(ct.sum(axis=1), axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.3, 1]})

    im = axes[0].imshow(ct_norm.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    axes[0].set_xticks(range(ct_norm.shape[1]))
    axes[0].set_xticklabels(ct_norm.columns, fontsize=8)
    axes[0].set_yticks(range(ct_norm.shape[0]))
    axes[0].set_yticklabels(ct_norm.index, fontsize=9)
    axes[0].set_xlabel("Consensus domain")
    axes[0].set_ylabel("Reference anatomical region")
    axes[0].set_title("A  Consensus domains recover known anatomy")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04).set_label("fraction of region")

    region_s3 = pd.DataFrame({"region": ref, "s3": s3["s3_score"].values})
    agg = region_s3.groupby("region")["s3"].agg(["mean", "std", "count"]).sort_values("mean")
    axes[1].barh(range(len(agg)), agg["mean"], xerr=agg["std"].fillna(0),
                 color="#1B9E77", ecolor="#555555", capsize=2)
    axes[1].set_yticks(range(len(agg)))
    axes[1].set_yticklabels([f"{r} (n={int(c)})" for r, c in zip(agg.index, agg["count"])],
                            fontsize=8)
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Mean S3")
    axes[1].set_title("B  Stability by anatomical region")

    fig.suptitle("Extended Data 7 | Validation against reference anatomy",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ==============================================================================
# Extended Data 8 -- cross-dataset generalization
# ==============================================================================
def ed8_cross_dataset(
    datasets: List[Dict[str, Any]],
    output_path: Optional[Path] = None,
) -> None:
    """
    Compare S3 across independent datasets to show the instability pattern
    generalizes.

    ``datasets`` is a list of {label, results_dir}. For each we load the per-spot
    S3 table and the canonical run's geometry. Top row: S3 stability maps; bottom
    row: overlaid S3 distributions and robust/fragile fractions per dataset.
    """
    import anndata as ad

    loaded = []
    for d in datasets:
        rd = Path(d["results_dir"])
        per_spot = pd.read_csv(rd / "aggregate" / "per_spot_stability.csv")
        matches = sorted((rd / "runs").glob("spec_00_*/adata_clustered.h5ad"))
        adata = ad.read_h5ad(matches[0]) if matches else None
        loaded.append({"label": d["label"], "per_spot": per_spot, "adata": adata})

    n = len(loaded)
    fig = plt.figure(figsize=(5.5 * n, 10))
    gs = gridspec.GridSpec(2, n, height_ratios=[1.6, 1], hspace=0.3, wspace=0.25)

    for j, d in enumerate(loaded):
        ax = fig.add_subplot(gs[0, j])
        if d["adata"] is not None:
            s3 = d["per_spot"].set_index("spot_id").reindex(d["adata"].obs_names)["s3_score"]
            sc = _spatial_scatter(ax, d["adata"].obsm["spatial"], s3.values, categorical=False,
                                  title=d["label"], cmap=S3_CMAP, vmin=0, vmax=1, s=8)
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label("S3")

    axd = fig.add_subplot(gs[1, :])
    for j, d in enumerate(loaded):
        vals = d["per_spot"]["s3_score"].values
        robust = 100 * np.mean(vals >= 0.8)
        median = np.median(vals)
        axd.hist(vals, bins=30, histtype="step", lw=2,
                 color=CLUSTER_COLORS[j],
                 label=f"{d['label']} (median {median:.2f}, {robust:.0f}% robust)")
    axd.axvline(0.8, color="#1A9850", ls="--", lw=1)
    axd.axvline(0.5, color="#D73027", ls="--", lw=1)
    axd.set_xlabel("Per-spot S3")
    axd.set_ylabel("Spots")
    axd.set_title("S3 distributions are conserved across datasets")
    axd.legend(fontsize=9)

    fig.suptitle("Extended Data 8 | Generalization across independent datasets",
                 fontsize=14, fontweight="bold")
    _save(fig, output_path)


# ==============================================================================
# Extended Data 9 -- interactive multiverse explorer (Plotly HTML)
# ==============================================================================
def generate_interactive_explorer(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    output_path: Path,
) -> None:
    """Self-contained Plotly HTML: slider across specs recolors the tissue live."""
    import plotly.graph_objects as go

    frames, steps = [], []
    for r in runs:
        a = r["adata"]
        codes = pd.Categorical(a.obs["leiden"]).codes
        frames.append(go.Frame(
            name=r["name"],
            data=[go.Scatter(
                x=a.obsm["spatial"][:, 0], y=a.obsm["spatial"][:, 1], mode="markers",
                marker=dict(size=5, color=codes, colorscale="Turbo"),
                text=[f"cluster {c}" for c in a.obs["leiden"]],
            )],
        ))
        steps.append(dict(method="animate", label=f"{r['idx']}:{r['name']}",
                          args=[[r["name"]], dict(mode="immediate",
                                                  frame=dict(duration=0, redraw=True))]))

    first = runs[0]["adata"]
    fig = go.Figure(
        data=[go.Scatter(
            x=first.obsm["spatial"][:, 0], y=first.obsm["spatial"][:, 1], mode="markers",
            marker=dict(size=5, color=pd.Categorical(first.obs["leiden"]).codes, colorscale="Turbo"),
        )],
        frames=frames,
    )
    fig.update_yaxes(autorange="reversed", scaleanchor="x", scaleratio=1)
    fig.update_layout(
        title="MOSAIC multiverse explorer — slide across specifications",
        sliders=[dict(active=0, steps=steps, currentvalue=dict(prefix="spec: "))],
        width=800, height=800, template="plotly_white",
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path), include_plotlyjs="cdn")
    logger.info("Wrote interactive explorer to %s", output_path)


# --- legacy helpers kept for the simple entry points ---------------------------
def _save(fig, output_path: Optional[Path]) -> None:
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def plot_specification_curve(metadata_list, metric="n_clusters_found",
                             output_path=None, title="Specification Curve"):
    """Sorted metric across specs (kept for backwards compatibility)."""
    values = np.array([m.get(metric, np.nan) for m in metadata_list], dtype=float)
    sorted_vals = np.sort(values)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(range(len(sorted_vals)), sorted_vals, "o-", color="steelblue", lw=2, ms=6)
    ax.set_xlabel("Specification (sorted)")
    ax.set_ylabel(metric)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    _save(fig, output_path)


def generate_all_figures(results_dir: Path, output_dir: Path, verbose: bool = True) -> None:
    """Generate the full figure set from a completed sweep."""
    results_dir = Path(results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    per_spot = pd.read_csv(results_dir / "aggregate" / "per_spot_stability.csv")
    per_feature = pd.read_csv(results_dir / "aggregate" / "per_feature_stability.csv")
    runs = load_runs(results_dir)
    if not runs:
        raise FileNotFoundError("No completed runs in results/runs/")
    canonical = runs[0]["adata"]

    figure1_spec_montage(runs, per_spot, output_dir / "fig1_spec_montage.png")
    figure2_spec_curve_matrix(runs, per_spot, output_dir / "fig2_spec_curve_matrix.png")
    figure3_stability_summary(canonical, per_spot, output_dir / "fig3_stability_summary.png")
    ed4_spec_space_table(runs, output_dir / "ed4_spec_space.png")
    ed5_per_feature_stability(per_feature, runs, output_dir / "ed5_svg_stability.png")
    ed6_headline_instability(runs, per_spot, output_dir / "ed6_headline_instability.png")
    ed7_reference_validation(runs, per_spot, output_dir / "ed7_reference_validation.png")
    generate_interactive_explorer(runs, per_spot, output_dir / "ed9_explorer.html")

    if verbose:
        logger.info("All figures written to %s", output_dir)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    generate_all_figures(project / "results", project / "results" / "figures", verbose=True)


if __name__ == "__main__":
    main()
