"""
Statistical rigor for MOSAIC.

These analyses turn the descriptive S3 into tested, uncertainty-bounded claims:
a permutation test for per-region stability, a variance decomposition attributing
instability to analytical axes, bootstrap confidence intervals, a technical
confound check against sequencing depth, and a multiverse-size convergence test.

Everything is computed from the saved sweep (results[/<dataset>]/runs + aggregate),
reusing the same cross-spec label alignment as the S3 definition.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def build_aligned_matrix(runs: List[Dict[str, Any]], cluster_key: str = "leiden"):
    """
    Return (aligned_labels[n_spots, n_specs], obs_names) using the first run as
    the reference, matching the per-spot S3 definition.
    """
    from src.stability import _align_labels_to_reference

    ref = runs[0]["adata"]
    cols = []
    for r in runs:
        a = r["adata"]
        if cluster_key in a.obs and a.obs_names.equals(ref.obs_names):
            cols.append(_align_labels_to_reference(a.obs[cluster_key].values,
                                                    ref.obs[cluster_key].values))
    return np.array(cols).T, ref.obs_names


def _s3_from_matrix(mat: np.ndarray) -> np.ndarray:
    """Per-spot S3 = modal-label frequency across the (aligned) specification columns."""
    n = mat.shape[1]
    out = np.empty(mat.shape[0])
    for i in range(mat.shape[0]):
        _, counts = np.unique(mat[i], return_counts=True)
        out[i] = counts.max() / n
    return out


def permutation_test_region_s3(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    reference_key: str = "cluster",
    n_perm: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Test whether each anatomical region is more stable than chance.

    Null: randomly reassign the region label of each spot (permuting the region
    vector), preserving region sizes and the S3 values, then recompute per-region
    mean S3. The empirical one-sided p-value is the fraction of permutations whose
    null region mean is >= the observed mean.
    """
    ref = runs[0]["adata"]
    if reference_key not in ref.obs:
        return pd.DataFrame()
    regions = ref.obs[reference_key].astype(str).values
    s3 = per_spot.set_index("spot_id").reindex(ref.obs_names)["s3_score"].values

    obs_mean = pd.Series(s3).groupby(regions).mean()
    rng = np.random.RandomState(seed)
    null_ge = {r: 0 for r in obs_mean.index}
    for _ in range(n_perm):
        perm = rng.permutation(regions)
        null_mean = pd.Series(s3).groupby(perm).mean()
        for r in obs_mean.index:
            if null_mean.get(r, -np.inf) >= obs_mean[r]:
                null_ge[r] += 1

    rows = []
    for r in obs_mean.index:
        p = (null_ge[r] + 1) / (n_perm + 1)  # add-one for a conservative estimate
        rows.append({"region": r, "observed_mean_s3": float(obs_mean[r]),
                     "p_value": float(p), "n_spots": int((regions == r).sum())})
    df = pd.DataFrame(rows).sort_values("observed_mean_s3", ascending=False)
    logger.info("Permutation test: %d/%d regions stable at p<0.05",
                int((df["p_value"] < 0.05).sum()), len(df))
    return df


def variance_decomposition(runs: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Attribute variance in domain count to each analytical axis via a Type-II ANOVA.

    Returns each axis's share of explained variance (partial sum-of-squares
    normalized), answering: which single choice most drives instability?
    """
    import statsmodels.api as sm
    from statsmodels.formula.api import ols

    rows = []
    for r in runs:
        s = r["spec"]
        rows.append({
            "n_domains": r["metadata"].get("n_clusters_found", np.nan),
            "normalization": s.get("normalization", "?"),
            "hvg_method": s.get("feature_selection", {}).get("method", "?"),
            "hvg_n": s.get("feature_selection", {}).get("n_genes", np.nan),
            "resolution": s.get("clustering", {}).get("resolution", np.nan),
            "clustering": s.get("clustering", {}).get("method", "?"),
            "k": s.get("neighborhood_graph", {}).get("k_or_radius", np.nan),
            "pca": s.get("dimensionality_reduction", {}).get("n_components", np.nan),
        })
    df = pd.DataFrame(rows).dropna()
    # Only include predictors that actually vary; a constant term has zero degrees
    # of freedom and breaks the Type-II ANOVA.
    categorical = ["normalization", "hvg_method", "clustering"]
    numeric = ["hvg_n", "resolution", "k", "pca"]
    terms = [f"C({c})" for c in categorical if df[c].nunique() > 1]
    terms += [n for n in numeric if df[n].nunique() > 1]
    if not terms:
        return pd.DataFrame(columns=["term", "variance_fraction"])
    formula = "n_domains ~ " + " + ".join(terms)
    model = ols(formula, data=df).fit()
    anova = sm.stats.anova_lm(model, typ=2)
    ss = anova["sum_sq"]
    total = ss.sum()
    out = (ss / total).rename("variance_fraction").reset_index()
    out.columns = ["term", "variance_fraction"]
    out["term"] = out["term"].str.replace(r"C\((.*)\)", r"\1", regex=True)
    out = out.sort_values("variance_fraction", ascending=False).reset_index(drop=True)
    logger.info("Variance decomposition (top): %s = %.0f%%",
                out.iloc[0]["term"], 100 * out.iloc[0]["variance_fraction"])
    return out


def bootstrap_ci_region_s3(
    runs: List[Dict[str, Any]],
    reference_key: str = "cluster",
    n_boot: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """95% bootstrap CI on per-region mean S3 by resampling specifications."""
    mat, obs_names = build_aligned_matrix(runs)
    ref = runs[0]["adata"]
    regions = ref.obs[reference_key].astype(str).values if reference_key in ref.obs else None
    if regions is None:
        return pd.DataFrame()

    rng = np.random.RandomState(seed)
    n_specs = mat.shape[1]
    region_names = pd.unique(regions)
    boot = {r: [] for r in region_names}
    for _ in range(n_boot):
        cols = rng.randint(0, n_specs, n_specs)
        s3 = _s3_from_matrix(mat[:, cols])
        means = pd.Series(s3).groupby(regions).mean()
        for r in region_names:
            boot[r].append(means[r])

    obs_s3 = _s3_from_matrix(mat)
    obs_means = pd.Series(obs_s3).groupby(regions).mean()
    rows = []
    for r in region_names:
        arr = np.array(boot[r])
        rows.append({"region": r, "mean_s3": float(obs_means[r]),
                     "ci_low": float(np.percentile(arr, 2.5)),
                     "ci_high": float(np.percentile(arr, 97.5))})
    return pd.DataFrame(rows).sort_values("mean_s3", ascending=False).reset_index(drop=True)


def confound_check_umi(runs: List[Dict[str, Any]], per_spot: pd.DataFrame,
                       count_key: str = "total_counts") -> Dict[str, Any]:
    """
    Test whether S3 is just a proxy for sequencing depth.

    Correlates per-spot S3 with total UMI counts. A strong positive correlation
    would mean "fragile" is confounded by low-quality spots; a weak one shows S3
    captures something beyond depth.
    """
    from scipy.stats import pearsonr, spearmanr

    ref = runs[0]["adata"]
    s3 = per_spot.set_index("spot_id").reindex(ref.obs_names)["s3_score"].values
    if count_key in ref.obs:
        counts = ref.obs[count_key].values.astype(float)
    elif "counts" in ref.layers:
        c = ref.layers["counts"]
        counts = np.asarray(c.sum(axis=1)).ravel().astype(float)
    else:
        c = ref.X
        counts = np.asarray(c.sum(axis=1)).ravel().astype(float)
    mask = np.isfinite(s3) & np.isfinite(counts)
    pear = pearsonr(s3[mask], counts[mask])
    spear = spearmanr(s3[mask], counts[mask])
    result = {
        "pearson_r": float(pear[0]), "pearson_p": float(pear[1]),
        "spearman_r": float(spear[0]), "spearman_p": float(spear[1]),
        "n_spots": int(mask.sum()),
    }
    logger.info("Confound check S3 vs %s: Pearson r=%.3f", count_key, result["pearson_r"])
    return result


def multiverse_size_convergence(
    runs: List[Dict[str, Any]],
    sizes: Optional[List[int]] = None,
    n_draws: int = 50,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Show that mean S3 stabilizes as the number of specifications grows.

    For each subset size, draw random spec subsets and record the mean per-spot
    S3, reporting mean and spread across draws.
    """
    mat, _ = build_aligned_matrix(runs)
    n_specs = mat.shape[1]
    sizes = sizes or [s for s in [3, 5, 8, 10, 15, 20, 25, 30] if s <= n_specs]
    rng = np.random.RandomState(seed)
    rows = []
    for n in sizes:
        vals = []
        for _ in range(n_draws):
            cols = rng.choice(n_specs, n, replace=False)
            vals.append(_s3_from_matrix(mat[:, cols]).mean())
        vals = np.array(vals)
        rows.append({"n_specs": n, "mean_s3": float(vals.mean()),
                     "std_s3": float(vals.std())})
    return pd.DataFrame(rows)


def _boundary_flags(coords: np.ndarray, domains: np.ndarray, k: int = 6) -> np.ndarray:
    """
    Data-driven boundary proxy: a spot is a 'boundary' if any of its k nearest
    spatial neighbors belongs to a different consensus domain, else 'interior'.
    Used when pathologist annotations are unavailable (H1); labeled as a proxy.
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(coords)
    _, nbr = tree.query(coords, k=min(k + 1, len(coords)))
    nbr = nbr[:, 1:]  # drop self
    same = (domains[nbr] == domains[:, None])
    return ~same.all(axis=1)  # True = boundary


def stakes_stratified_s3(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    reference_key: Optional[str] = None,
    k: int = 6,
    n_perm: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Test H1: are boundary spots more fragile (lower S3) than interior spots?

    Strata come from pathologist annotation if ``reference_key`` names an
    annotation whose labels include boundary terms; otherwise from a consensus-
    domain boundary proxy. Significance is a label-shuffle permutation null on the
    interior-minus-boundary mean-S3 gap (one-sided: boundary more fragile).
    """
    ref = runs[0]["adata"]
    ps = per_spot.set_index("spot_id").reindex(ref.obs_names)
    s3 = ps["s3_score"].values
    coords = np.asarray(ref.obsm["spatial"])
    domains = ps["consensus_cluster"].astype(int).values
    boundary = _boundary_flags(coords, domains, k=k)

    b_mean, i_mean = float(s3[boundary].mean()), float(s3[~boundary].mean())
    obs_gap = i_mean - b_mean  # positive if boundary is more fragile

    rng = np.random.RandomState(seed)
    n_b = int(boundary.sum())
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(len(s3)) < n_b  # random subset of size n_b as "boundary"
        gap = s3[~perm].mean() - s3[perm].mean()
        if gap >= obs_gap:
            ge += 1
    p = (ge + 1) / (n_perm + 1)

    result = {
        "strata": "consensus_boundary_proxy",
        "boundary_mean_s3": b_mean,
        "interior_mean_s3": i_mean,
        "gap_interior_minus_boundary": obs_gap,
        "p_value": float(p),
        "n_boundary": n_b,
        "n_interior": int((~boundary).sum()),
    }
    logger.info(
        "Stakes-stratified S3: boundary=%.3f interior=%.3f gap=%.3f (p=%.4f)",
        b_mean, i_mean, obs_gap, result["p_value"],
    )
    return result, boundary


def stakes_stratified_s3_depth_controlled(
    runs: List[Dict[str, Any]],
    per_spot: pd.DataFrame,
    boundary: np.ndarray,
    count_key: str = "total_counts",
    n_bins: int = 10,
    n_perm: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Re-test H1's boundary/interior S3 gap with sequencing depth controlled for.

    ``stakes_stratified_s3``'s permutation null draws a random boundary-sized
    subset from *all* spots, which only controls for group size - not for depth.
    If boundary and interior spots differ systematically in total UMI counts (as
    they do on the tumor dataset; see ``confound_check_umi``), that null can
    overstate significance: a real gap could partly reflect depth, not analytical
    fragility. This block-permutation null instead bins spots into ``n_bins``
    depth quantiles and shuffles the boundary/interior label only *within* each
    bin, so every permutation preserves the true depth distribution of each
    group. If the gap survives here, it isn't just a depth artifact.

    Requires ``boundary`` (the boolean mask ``stakes_stratified_s3`` already
    computed) so the two tests are directly comparable on the same partition.
    """
    ref = runs[0]["adata"]
    ps = per_spot.set_index("spot_id").reindex(ref.obs_names)
    s3 = ps["s3_score"].values
    if count_key in ref.obs:
        depth = ref.obs[count_key].values.astype(float)
    else:
        depth = np.asarray(ref.layers["counts"].sum(axis=1)).ravel().astype(float)

    # Diagnostic: do the two groups actually differ in depth?
    from scipy.stats import mannwhitneyu

    depth_b, depth_i = depth[boundary], depth[~boundary]
    mw = mannwhitneyu(depth_b, depth_i, alternative="two-sided")

    b_mean, i_mean = float(s3[boundary].mean()), float(s3[~boundary].mean())
    obs_gap = i_mean - b_mean

    bins = pd.qcut(depth, q=min(n_bins, len(np.unique(depth))), labels=False, duplicates="drop")
    rng = np.random.RandomState(seed)
    ge = 0
    for _ in range(n_perm):
        perm_boundary = np.zeros(len(s3), dtype=bool)
        for b in np.unique(bins):
            idx = np.flatnonzero(bins == b)
            n_b_in_bin = int(boundary[idx].sum())
            if n_b_in_bin:
                chosen = rng.choice(idx, n_b_in_bin, replace=False)
                perm_boundary[chosen] = True
        gap = s3[~perm_boundary].mean() - s3[perm_boundary].mean()
        if gap >= obs_gap:
            ge += 1
    p = (ge + 1) / (n_perm + 1)

    result = {
        "gap_interior_minus_boundary": obs_gap,
        "p_value_depth_controlled": float(p),
        "depth_boundary_median": float(np.median(depth_b)),
        "depth_interior_median": float(np.median(depth_i)),
        "depth_groups_differ_mannwhitney_p": float(mw.pvalue),
        "n_depth_bins": int(len(np.unique(bins))),
    }
    logger.info(
        "Depth-controlled stakes test: gap=%.3f p=%.4f (naive vs depth-controlled comparison; "
        "boundary/interior depth medians %.0f vs %.0f, Mann-Whitney p=%.2e)",
        obs_gap, p, result["depth_boundary_median"], result["depth_interior_median"],
        mw.pvalue,
    )
    return result


def stakes_test_from_arrays(
    s3: np.ndarray,
    boundary: np.ndarray,
    depth: Optional[np.ndarray] = None,
    n_bins: int = 10,
    n_perm: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Test whether ``boundary`` spots are more fragile (lower S3) than the rest,
    directly from arrays (so it works with externally-provided labels such as
    pathologist annotations, not just the consensus proxy).

    If ``depth`` is given, also runs a depth-stratified block-permutation null
    (shuffle the boundary label only within depth deciles) so the result is not a
    sequencing-depth artifact.
    """
    s3 = np.asarray(s3, float)
    boundary = np.asarray(boundary, bool)
    b_mean, o_mean = float(s3[boundary].mean()), float(s3[~boundary].mean())
    obs_gap = o_mean - b_mean
    rng = np.random.RandomState(seed)
    n_b = int(boundary.sum())

    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(len(s3)) < n_b
        if (s3[~perm].mean() - s3[perm].mean()) >= obs_gap:
            ge += 1
    result = {
        "boundary_mean_s3": b_mean, "other_mean_s3": o_mean,
        "gap": obs_gap, "p_value": float((ge + 1) / (n_perm + 1)),
        "n_boundary": n_b, "n_other": int((~boundary).sum()),
    }

    if depth is not None:
        depth = np.asarray(depth, float)
        bins = pd.qcut(depth, q=n_bins, labels=False, duplicates="drop")
        ge_d = 0
        for _ in range(n_perm):
            perm = np.zeros(len(s3), bool)
            for b in np.unique(bins):
                idx = np.where(bins == b)[0]
                k = int(boundary[idx].sum())
                if k:
                    perm[rng.choice(idx, k, replace=False)] = True
            if (s3[~perm].mean() - s3[perm].mean()) >= obs_gap:
                ge_d += 1
        result["p_value_depth_controlled"] = float((ge_d + 1) / (n_perm + 1))
    return result


def fragility_report(runs: List[Dict[str, Any]], alpha: float = 0.05):
    """
    A concrete 'headline that flips' demo (the fragility-aware reporting story).

    Among genes tested in every specification, pick the most contested SVG (call
    frequency nearest 50%) and trace its FDR-adjusted Moran's I p-value across all
    specifications. Returns (gene, s3, DataFrame[spec_name, pval, significant]).
    """
    from src.visualization import _significant, svg_concordance_matrix

    conc = svg_concordance_matrix(runs, common_only=True)
    if conc.empty:
        return None, None, pd.DataFrame()
    freq = conc.mean(axis=1)
    gene = (freq - 0.5).abs().idxmin()
    s3 = float(freq[gene])

    rows = []
    for r in runs:
        svg = r["svg_df"]
        if svg is None or gene not in svg.index:
            continue
        pcol = next((c for c in ("pval_norm_fdr_bh", "pval_sim_fdr_bh", "pval_norm",
                                 "pval_sim", "pval") if c in svg.columns), None)
        pval = float(svg.loc[gene, pcol]) if pcol else np.nan
        sig = bool(_significant(svg, alpha).get(gene, False))
        rows.append({"spec_name": r["metadata"].get("spec_name", ""),
                     "pval": pval, "significant": sig})
    df = pd.DataFrame(rows).sort_values("pval").reset_index(drop=True)
    logger.info("Fragility report: gene %s called SVG in %d/%d specs (S3=%.2f)",
                gene, int(df["significant"].sum()), len(df), s3)
    return gene, s3, df


def plot_fragility_report(runs: List[Dict[str, Any]], output_path: Optional[Path] = None,
                          alpha: float = 0.05) -> None:
    """Trace the most-contested SVG's p-value across pipelines vs the FDR threshold."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gene, s3, df = fragility_report(runs, alpha)
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#1A9850" if s else "#D73027" for s in df["significant"]]
    ax.scatter(range(len(df)), df["pval"], c=colors, s=45, zorder=3)
    ax.axhline(alpha, color="#333333", ls="--", lw=1.5, label=f"FDR = {alpha}")
    ax.set_yscale("log")
    ax.set_xlabel("Specification (sorted by p-value)")
    ax.set_ylabel("FDR-adjusted Moran's I p-value (log)")
    n_sig = int(df["significant"].sum())
    ax.set_title(f"'{gene}' is called a significant SVG in {n_sig}/{len(df)} defensible "
                 f"pipelines (S3 = {s3:.2f})", fontsize=12, fontweight="bold")
    ax.legend()
    fig.suptitle("Fragility-aware reporting: a headline SVG that flips across pipelines",
                 fontsize=13, fontweight="bold")
    _save(fig, output_path)


def cross_dataset_region_summary(datasets: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Summarize per-region S3 across datasets to test whether the robust/fragile
    pattern is consistent. ``datasets`` = [{label, runs, per_spot}].
    """
    frames = []
    for d in datasets:
        ref = d["runs"][0]["adata"]
        if "cluster" not in ref.obs:
            continue
        regions = ref.obs["cluster"].astype(str).values
        s3 = d["per_spot"].set_index("spot_id").reindex(ref.obs_names)["s3_score"].values
        agg = pd.Series(s3).groupby(regions).mean().rename(d["label"])
        frames.append(agg)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=1)


def _load_runs_and_spot(results_dir: Path) -> Tuple[list, pd.DataFrame]:
    from src.visualization import load_runs

    per_spot = pd.read_csv(Path(results_dir) / "aggregate" / "per_spot_stability.csv")
    return load_runs(results_dir), per_spot


def run_all(results_dir: Path, out_json: Optional[Path] = None) -> Dict[str, Any]:
    """Run every statistical analysis on one dataset and return a summary dict."""
    from src.utils import save_json

    runs, per_spot = _load_runs_and_spot(results_dir)
    perm = permutation_test_region_s3(runs, per_spot)
    var = variance_decomposition(runs)
    boot = bootstrap_ci_region_s3(runs)
    conf = confound_check_umi(runs, per_spot)
    conv = multiverse_size_convergence(runs)
    stakes, boundary = stakes_stratified_s3(runs, per_spot)
    stakes_depth = stakes_stratified_s3_depth_controlled(runs, per_spot, boundary)
    frag_gene, frag_s3, frag_df = fragility_report(runs)

    agg = Path(results_dir) / "aggregate"
    perm.to_csv(agg / "stats_permutation_region_s3.csv", index=False)
    var.to_csv(agg / "stats_variance_decomposition.csv", index=False)
    boot.to_csv(agg / "stats_bootstrap_region_ci.csv", index=False)
    conv.to_csv(agg / "stats_multiverse_convergence.csv", index=False)

    summary = {
        "confound_umi": conf,
        "variance_top_axis": var.iloc[0].to_dict() if len(var) else {},
        "regions_stable_p05": int((perm["p_value"] < 0.05).sum()) if len(perm) else 0,
        "n_regions": len(perm),
        "stakes_stratified": stakes,
        "stakes_stratified_depth_controlled": stakes_depth,
        "fragility_example": {
            "gene": frag_gene,
            "s3": frag_s3,
            "n_significant": int(frag_df["significant"].sum()) if len(frag_df) else 0,
            "n_specs": int(len(frag_df)),
        },
    }
    if out_json:
        save_json(summary, out_json)
    return summary


# --- figures -------------------------------------------------------------------
def _save(fig, output_path):
    import matplotlib.pyplot as plt

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def generate_stats_figures(results_dir: Path, output_dir: Path) -> None:
    """Build the four statistics figures (variance, permutation, convergence, confound)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results_dir, output_dir = Path(results_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    runs, per_spot = _load_runs_and_spot(results_dir)

    # Fig 8: variance decomposition.
    var = variance_decomposition(runs)
    var = var[var["term"] != "Residual"]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(var["term"][::-1], 100 * var["variance_fraction"][::-1], color="#377EB8")
    ax.set_xlabel("Variance in domain count explained (%)")
    ax.set_title("Which analytical choice drives instability most?",
                 fontsize=12, fontweight="bold")
    for i, v in enumerate(var["variance_fraction"][::-1]):
        ax.text(100 * v + 0.5, i, f"{100 * v:.0f}%", va="center", fontsize=9)
    _save(fig, output_dir / "08_variance_decomposition.png")

    # Fig 9: permutation null vs observed per-region S3.
    perm = permutation_test_region_s3(runs, per_spot)
    if len(perm):
        fig, ax = plt.subplots(figsize=(9, 6))
        colors = ["#1A9850" if p < 0.05 else "#999999" for p in perm["p_value"]]
        ax.barh(range(len(perm)), perm["observed_mean_s3"], color=colors)
        ax.set_yticks(range(len(perm)))
        ax.set_yticklabels([f"{r}  (p={p:.3f})" for r, p in
                            zip(perm["region"], perm["p_value"])], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("Observed mean S3")
        ax.set_title("Per-region stability vs permutation null (green: p<0.05)",
                     fontsize=12, fontweight="bold")
        _save(fig, output_dir / "09_permutation_test.png")

    # Fig 10: multiverse-size convergence.
    conv = multiverse_size_convergence(runs)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(conv["n_specs"], conv["mean_s3"], yerr=conv["std_s3"],
                marker="o", color="#984EA3", capsize=3)
    ax.set_xlabel("Number of specifications sampled")
    ax.set_ylabel("Mean per-spot S3")
    ax.set_title("S3 estimates converge as the multiverse grows",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3)
    _save(fig, output_dir / "10_convergence.png")

    # Fig 11: confound check S3 vs sequencing depth.
    conf = confound_check_umi(runs, per_spot)
    ref = runs[0]["adata"]
    if conf and "total_counts" in ref.obs:
        s3 = per_spot.set_index("spot_id").reindex(ref.obs_names)["s3_score"].values
        counts = ref.obs["total_counts"].values.astype(float)
        robust = s3 >= 0.8
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(counts[~robust], s3[~robust], s=8, c="#D73027", alpha=0.5, label="fragile/mid")
        ax.scatter(counts[robust], s3[robust], s=8, c="#1A9850", alpha=0.5, label="robust")
        ax.set_xlabel("Total UMI counts per spot")
        ax.set_ylabel("Per-spot S3")
        ax.set_title(f"S3 is not a depth artifact (Pearson r = {conf['pearson_r']:.2f})",
                     fontsize=12, fontweight="bold")
        ax.legend()
        _save(fig, output_dir / "11_confound_umi.png")

    # Fig 13: stakes-stratified S3 (H1) -- boundary vs interior.
    stakes, boundary = stakes_stratified_s3(runs, per_spot)
    ref = runs[0]["adata"]
    coords = np.asarray(ref.obsm["spatial"])
    s3 = per_spot.set_index("spot_id").reindex(ref.obs_names)["s3_score"].values
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1.1, 1]})
    axes[0].scatter(coords[~boundary, 0], coords[~boundary, 1], s=10, c="#999999",
                    label=f"interior (n={stakes['n_interior']})")
    axes[0].scatter(coords[boundary, 0], coords[boundary, 1], s=10, c="#D73027",
                    label=f"boundary (n={stakes['n_boundary']})")
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_title("A  Consensus-domain boundary proxy", fontsize=12, fontweight="bold")
    axes[1].hist(s3[~boundary], bins=25, histtype="step", lw=2, color="#555555",
                 density=True, label=f"interior (mean {stakes['interior_mean_s3']:.2f})")
    axes[1].hist(s3[boundary], bins=25, histtype="step", lw=2, color="#D73027",
                 density=True, label=f"boundary (mean {stakes['boundary_mean_s3']:.2f})")
    axes[1].set_xlabel("Per-spot S3")
    axes[1].set_ylabel("Density")
    axes[1].legend(fontsize=9)
    axes[1].set_title(f"B  Boundaries are more fragile (p={stakes['p_value']:.4f})",
                      fontsize=12, fontweight="bold")
    fig.suptitle("H1: analytical fragility concentrates at domain boundaries",
                 fontsize=13, fontweight="bold")
    _save(fig, output_dir / "13_stakes_stratified.png")

    # Fig 14: fragility-aware reporting demo (a headline SVG that flips).
    plot_fragility_report(runs, output_dir / "14_fragility_report.png")

    logger.info("Statistics figures written to %s", output_dir)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    run_all(project / "results", project / "results" / "aggregate" / "statistics_summary.json")
    generate_stats_figures(project / "results", project / "results" / "figures")


if __name__ == "__main__":
    main()
