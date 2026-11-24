"""
Populate paper_ieee/numbers.tex from the aggregate sweep results.

Run after scripts/main.py:  PYTHONPATH=. python scripts/report_numbers.py
Reads results/aggregate/*, writes LaTeX \\newcommand macros used by paper_ieee/main.tex.
"""

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent
AGG = PROJECT / "results" / "aggregate"
OUT = PROJECT / "paper_ieee" / "numbers.tex"


def _fmt(x, pct=False, digits=1):
    if pct:
        return f"{100 * x:.0f}\\%"
    return f"{x:.{digits}f}"



def _pct(x):
    return f"{100 * x:.0f}\\%"


def _load(p):
    p = PROJECT / p
    return json.loads(p.read_text()) if p.exists() else None


def _extra_macros() -> dict:
    """Macros that used to be added to numbers.tex by hand."""
    import importlib.util

    import numpy as np

    m = {}
    # Pathology gaps (scripts/pathology_gaps.py).
    pg = _load("results/breast_cancer/aggregate/pathology_margin_vs_invasive.json")
    if pg:
        m["PathMarginInvasiveGap"] = f"{pg['margin_vs_invasive_gap']:.2f}"
        m["PathMarginInvasiveCILow"] = f"{pg['margin_vs_invasive_ci95'][0]:.3f}"
        m["PathMarginInvasiveCIHigh"] = f"{pg['margin_vs_invasive_ci95'][1]:.3f}"
        m["PathMarginRestGap"] = f"{pg['margin_vs_rest_gap']:.2f}"
        m["PathMarginRestCILow"] = f"{pg['margin_vs_rest_ci95'][0]:.3f}"
        m["PathMarginRestCIHigh"] = f"{pg['margin_vs_rest_ci95'][1]:.3f}"
    # Gene-level selection frequency vs call consistency.
    pf = PROJECT / "results" / "aggregate" / "per_feature_stability.csv"
    if pf.exists():
        s3 = pd.read_csv(pf)["s3_score"]
        m["GeneSelFreqLowFrac"] = _pct((s3 < 0.5).mean())
        m["GeneConsistentFrac"] = _pct((np.maximum(s3, 1 - s3) > 0.8).mean())
    # Single-pipeline predictor (mosaic.predict_fragility).
    fr = _load("results/aggregate/fragility_prediction.json")
    fs = _load("results/aggregate/fragility_prediction_spatial_cv.json")
    fb = _load("results/aggregate/fragility_prediction_buffered_spatial_cv.json")
    if fr:
        m["FragPredAUCRandomCILow"] = f"{fr['significance']['ci95_low']:.3f}"
        m["FragPredAUCRandomCIHigh"] = f"{fr['significance']['ci95_high']:.3f}"
        m["FragPredPermP"] = "$p\\le0.001$"
        m["FragPredBestFeatAUC"] = f"{max(fr['per_feature_auc'].values()):.2f}"
        m["FragPredDepthAUC"] = f"{fr['per_feature_auc']['depth']:.2f}"
        m["FragPredNFragile"] = str(fr["n_fragile"])
    if fs:
        m["FragPredAUCSpatial"] = f"{fs['auc_combined_spatial_cv']:.2f}"
        m["FragPredAUCSpatialCILow"] = f"{fs['significance']['ci95_low']:.3f}"
        m["FragPredAUCSpatialCIHigh"] = f"{fs['significance']['ci95_high']:.3f}"
        m["FragPredSpatialBlocks"] = str(fs["n_blocks"])
    if fb:
        m["FragPredAUCBuffered"] = f"{fb['auc_combined_buffered_spatial_cv']:.2f}"
        m["FragPredAUCBufferedCILow"] = f"{fb['significance']['ci95_low']:.3f}"
        m["FragPredAUCBufferedCIHigh"] = f"{fb['significance']['ci95_high']:.3f}"
        m["FragPredBufferMultiple"] = f"{fb['buffer_multiple']:.0f}"
        m["FragPredBufferDroppedPerFold"] = f"{fb['mean_train_spots_dropped_per_fold']:.0f}"
    # Reference dependence of per-spot S3 (mosaic.reference_robustness).
    rr = _load("results/aggregate/reference_robustness.json")
    if rr:
        m["RefRobustN"] = str(rr["n_references"])
        m["RefRobustMeanPearson"] = f"{rr['mean_pairwise_pearson_r']:.2f}"
        m["RefRobustMinPearson"] = f"{rr['min_pairwise_pearson_r']:.2f}"
        m["RefRobustMeanPerSpotStd"] = f"{rr['mean_per_spot_std_across_refs']:.3f}"
        m["RefRobustAggMeanStd"] = f"{rr['s3_across_refs_std_of_means']:.3f}"
        lo, hi = rr.get("s3_across_refs_min"), rr.get("s3_across_refs_max")
        if lo is not None and hi is not None:
            m["RefRobustAggMeanRange"] = f"{lo:.2f}, {hi:.2f}"
    # Reference-free consensus selection (mosaic.consensus_selection).
    cs = _load("results/aggregate/consensus_selection_blind.json")
    if cs:
        m["ConsensusBlindLinkage"] = cs["blind_linkage"]
        m["ConsensusBlindK"] = str(cs["blind_k"])
        m["ConsensusBlindARI"] = f"{cs['blind_ari']:.3f}"
        m["ConsensusOracleARI"] = f"{cs['oracle_ari']:.3f}"
        m["ConsensusBlindOracleGap"] = f"{cs['ari_gap_blind_vs_oracle']:.3f}"
        m["ConsensusSilARISpearman"] = f"{cs['silhouette_ari_spearman']:.2f}"
        m["ConsensusNSettings"] = str(cs["n_settings_evaluated"])
    cr = _load("results/aggregate/consensus_recovery_summary.json")
    if cr and "uncertainty" in cr:
        m["ConsensusSthreeCorrect"] = f"{cr['uncertainty']['mean_s3_correct']:.2f}"
        m["ConsensusSthreeWrong"] = f"{cr['uncertainty']['mean_s3_incorrect']:.2f}"
    # Full factorial (results/full_*) and categorical variance shares.
    spec = importlib.util.spec_from_file_location(
        "variance_categorical", PROJECT / "scripts" / "variance_categorical.py")
    vc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vc)

    def cat(prefix, path):
        if not (PROJECT / path).exists():
            return
        base, ci, n = vc.shares_with_ci(PROJECT / path)
        for term, name in (("res", "Res"), ("k", "K"), ("norm", "Norm"),
                           ("hvgm", "Hvg"), ("Residual", "Resid")):
            if term in base:
                m[f"{prefix}{name}"] = _pct(base[term])
                lo, hi = ci[term]
                m[f"{prefix}{name}CI"] = f"{100 * lo:.0f}--{100 * hi:.0f}\\%"
        return n

    cat("CatHne", "results/aggregate/specification_curves.json")
    cat("CatBreast", "results/breast_cancer/aggregate/specification_curves.json")
    n_full = cat("FullHne", "results/full_hne/aggregate/specification_curves.json")
    cat("FullBreast", "results/full_breast_cancer/aggregate/specification_curves.json")
    if n_full:
        m["FullNspecs"] = str(n_full)
    for tag, d in (("Hne", "full_hne"), ("Breast", "full_breast_cancer")):
        pa = _load(f"results/{d}/aggregate/partition_concordance_ari.json")
        st = _load(f"results/{d}/aggregate/statistics_summary.json")
        ps = PROJECT / "results" / d / "aggregate" / "per_spot_stability.csv"
        if pa:
            m[f"FullAri{tag}"] = f"{pa['mean']:.2f}"
        if ps.exists():
            m[f"FullFragile{tag}"] = _pct((pd.read_csv(ps)["s3_score"] < 0.5).mean())
        if st:
            m[f"FullGap{tag}"] = f"{st['stakes_stratified']['gap_interior_minus_boundary']:.2f}"
    return m

def main() -> None:
    per_spot = pd.read_csv(AGG / "per_spot_stability.csv")
    per_feature = pd.read_csv(AGG / "per_feature_stability.csv")
    headline = json.loads((AGG / "headline_instability_n_clusters.json").read_text())
    specs = json.loads((AGG / "specification_curves.json").read_text())["metadata"]
    ari_path = AGG / "partition_concordance_ari.json"
    ari = json.loads(ari_path.read_text()) if ari_path.exists() else {}
    ref_path = AGG / "reference_concordance_ari.json"
    ref = json.loads(ref_path.read_text()) if ref_path.exists() else {}
    stats_path = AGG / "statistics_summary.json"
    stats = json.loads(stats_path.read_text()) if stats_path.exists() else {}
    down_path = AGG / "downstream_summary.json"
    down = json.loads(down_path.read_text()) if down_path.exists() else {}

    n_specs = len(specs)
    n_spots = len(per_spot)

    macros = {
        "Nspecs": str(n_specs),
        "Nspots": str(n_spots),
        "Ngenes": "18078",
        "ClustMin": f"{headline['min']:.0f}",
        "ClustMax": f"{headline['max']:.0f}",
        "ClustMean": _fmt(headline["mean"]),
        "ClustStd": _fmt(headline["std"]),
        "SpotSthreeMean": _fmt(per_spot["s3_score"].mean(), digits=2),
        "SpotStableFrac": _fmt((per_spot["s3_score"] >= 0.8).mean(), pct=True),
        "SpotFragileFrac": _fmt((per_spot["s3_score"] < 0.5).mean(), pct=True),
        "GeneRobustFrac": _fmt((per_feature["s3_score"] > 0.8).mean(), pct=True),
        "GeneFragileFrac": _fmt((per_feature["s3_score"] < 0.5).mean(), pct=True),
        "AriMean": _fmt(ari["mean"], digits=2) if ari else "TBD",
        "AriMin": _fmt(ari["min"], digits=2) if ari else "TBD",
        "AriMax": _fmt(ari["max"], digits=2) if ari else "TBD",
        "RefAriMean": _fmt(ref["mean"], digits=2) if ref else "TBD",
        "RefAriMax": _fmt(ref["max"], digits=2) if ref else "TBD",
    }
    # Per-axis variance shares (Type-II ANOVA) for the axes the paper discusses.
    for prefix, csv_path in (("Var", PROJECT / "results" / "aggregate" / "stats_variance_decomposition.csv"),
                             ("TumorVar", PROJECT / "results" / "breast_cancer" / "aggregate" / "stats_variance_decomposition.csv")):
        if csv_path.exists():
            vd = pd.read_csv(csv_path).set_index("term")["variance_fraction"]
            for term, name in (("k", "K"), ("normalization", "Norm"), ("Residual", "Resid")):
                if term in vd:
                    macros[f"{prefix}{name}Frac"] = _fmt(vd[term], pct=True)
    if stats:
        top = stats.get("variance_top_axis", {})
        conf = stats.get("confound_umi", {})
        macros["VarTopAxis"] = str(top.get("term", "TBD")).replace("_", " ")
        macros["VarTopFrac"] = _fmt(top.get("variance_fraction", 0), pct=True) if top else "TBD"
        macros["ConfoundR"] = _fmt(conf.get("pearson_r", 0), digits=2) if conf else "TBD"
        macros["RegionsStable"] = str(stats.get("regions_stable_p05", "TBD"))
        macros["NregionsTotal"] = str(stats.get("n_regions", "TBD"))
        st = stats.get("stakes_stratified", {})
        if st:
            macros["BoundaryS"] = _fmt(st.get("boundary_mean_s3", 0), digits=2)
            macros["InteriorS"] = _fmt(st.get("interior_mean_s3", 0), digits=2)
            macros["StakesGap"] = _fmt(st.get("gap_interior_minus_boundary", 0), digits=2)
            p = st.get("p_value", 1)
            macros["StakesP"] = "\\le 0.001" if p < 0.001 else _fmt(p, digits=3)
    if down:
        macros["DownstreamRangePP"] = _fmt(down.get("accuracy_range_pp", 0))
        macros["DownstreamAccMin"] = _fmt(down.get("accuracy_min", 0), digits=2)
        macros["DownstreamAccMax"] = _fmt(down.get("accuracy_max", 0), digits=2)
        if "brier_min" in down:
            macros["BrierMin"] = _fmt(down["brier_min"], digits=3)
            macros["BrierMax"] = _fmt(down["brier_max"], digits=3)

    # Multiverse consensus as a better estimator (anatomy recovery vs single pipelines).
    cons_path = AGG / "consensus_recovery_summary.json"
    if cons_path.exists():
        c = json.loads(cons_path.read_text())
        macros["ConsensusARI"] = _fmt(c.get("consensus_at_ref_k") or c.get("consensus_best_ari", 0), digits=2)
        macros["SingleMedianARI"] = _fmt(c.get("single_median", 0), digits=2)
        macros["SingleBestARI"] = _fmt(c.get("single_max", 0), digits=2)
        macros["SingleWorstARI"] = _fmt(c.get("single_min", 0), digits=2)
        macros["ConsensusBeatsMedianPct"] = _fmt(c.get("beats_median_frac", 0), pct=True)
        unc = c.get("uncertainty", {})
        if unc:
            macros["ConsensusUncAUC"] = _fmt(unc.get("auc_s3_predicts_correct", 0), digits=2)
            macros["ConsensusAcc"] = _fmt(unc.get("consensus_accuracy", 0), pct=True)

    # Single-pipeline fragility early-warning (AUC).
    fp = AGG / "fragility_prediction.json"
    if fp.exists():
        macros["FragPredAUC"] = _fmt(json.loads(fp.read_text()).get("auc_combined", 0), digits=2)

    # Second dataset: human breast-cancer tumor (H1/H2 on higher-stakes tissue).
    tumor_path = PROJECT / "results" / "breast_cancer" / "aggregate" / "statistics_summary.json"
    if tumor_path.exists():
        t = json.loads(tumor_path.read_text())
        ts = t.get("stakes_stratified", {})
        td = t.get("stakes_stratified_depth_controlled", {})
        tv = t.get("variance_top_axis", {})
        tc = t.get("confound_umi", {})
        if ts:
            macros["TumorBoundaryS"] = _fmt(ts.get("boundary_mean_s3", 0), digits=2)
            macros["TumorInteriorS"] = _fmt(ts.get("interior_mean_s3", 0), digits=2)
            macros["TumorGap"] = _fmt(ts.get("gap_interior_minus_boundary", 0), digits=2)
            tp = ts.get("p_value", 1)
            macros["TumorStakesP"] = "\\le 0.001" if tp < 0.001 else _fmt(tp, digits=3)
        if tc:
            macros["TumorConfoundR"] = _fmt(tc.get("pearson_r", 0), digits=2)
        if tv:
            macros["TumorVarTopAxis"] = str(tv.get("term", "")).replace("_", " ")
            macros["TumorVarTopFrac"] = _fmt(tv.get("variance_fraction", 0), pct=True)
        if td:
            dp = td.get("p_value_depth_controlled", 1)
            macros["TumorDepthP"] = "\\le 0.001" if dp < 0.001 else _fmt(dp, digits=3)

    # Independent pathologist-annotation validation of tumor H1 (SEDR labels).
    path_p = PROJECT / "results" / "breast_cancer" / "aggregate" / "pathology_stakes.json"
    if path_p.exists():
        pj = json.loads(path_p.read_text())
        cats = pj.get("per_category_mean_s3", {})
        macros["PathMarginS"] = _fmt(pj.get("margin_mean_s3", 0), digits=2)
        macros["PathInvasiveS"] = _fmt(cats.get("Invasive", 0), digits=2)
        macros["PathHealthyS"] = _fmt(cats.get("Healthy", 0), digits=2)
        pp = pj.get("p_value", 1)
        pdp = pj.get("p_value_depth_controlled", 1)
        macros["PathP"] = "\\le 0.001" if pp < 0.001 else _fmt(pp, digits=3)
        macros["PathDepthP"] = "\\le 0.001" if pdp < 0.001 else _fmt(pdp, digits=3)

    macros.update(_extra_macros())

    lines = [
        "% Auto-generated by scripts/report_numbers.py -- do not edit by hand.",
    ]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()]
    OUT.write_text("\n".join(lines) + "\n")

    print(f"Wrote {OUT} with {len(macros)} macros:")
    for k, v in macros.items():
        print(f"  {k} = {v}")


if __name__ == "__main__":
    if not AGG.exists():
        print("No aggregate results found; run scripts/main.py first.", file=sys.stderr)
        sys.exit(1)
    main()
