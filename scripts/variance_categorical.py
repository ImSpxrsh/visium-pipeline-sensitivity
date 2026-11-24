"""Type-II ANOVA of domain count with every axis categorical, plus 95% bootstrap
intervals from resampling pipelines. Usage: python scripts/variance_categorical.py SPEC_CURVES.json ...
"""
import json,sys,numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.formula.api import ols
def table(path):
    rows=[]
    for m in json.load(open(path))["metadata"]:
        s=m["spec"]
        rows.append(dict(y=m["n_clusters_found"],norm=s["normalization"],hvgm=s["feature_selection"]["method"],
          hvgn=s["feature_selection"]["n_genes"],res=s["clustering"]["resolution"],cl=s["clustering"]["method"],
          k=s["neighborhood_graph"]["k_or_radius"],pc=s["dimensionality_reduction"]["n_components"]))
    return pd.DataFrame(rows)
def shares(df):
    terms=[f"C({c})" for c in ["norm","hvgm","hvgn","res","cl","k","pc"] if df[c].nunique()>1]
    a=sm.stats.anova_lm(ols("y ~ "+" + ".join(terms),data=df).fit(),typ=2)["sum_sq"]
    sh=a/a.sum(); sh.index=[i.replace("C(","").replace(")","") for i in sh.index]; return sh
def shares_with_ci(path, n_boot=1000, seed=0):
    """Point shares and 95% bootstrap intervals (resampling pipelines)."""
    df = table(path)
    base = shares(df)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        try:
            boots.append(shares(df.sample(len(df), replace=True,
                                          random_state=int(rng.integers(1e9)))))
        except Exception:  # a resample can drop a factor level entirely
            pass
    b = pd.DataFrame(boots)
    ci = {t: (float(np.nanpercentile(b[t], 2.5)), float(np.nanpercentile(b[t], 97.5)))
          for t in base.index if t in b}
    return base, ci, len(df)


LABELS = {"res": "Resolution", "k": "Neighbors (k)", "norm": "Normalization",
          "hvgm": "HVG method", "hvgn": "HVG count", "pc": "PCs",
          "cl": "Cluster method", "Residual": "Residual"}


def plot(path, out, title=None):
    """Bar chart of categorical shares with 95% bootstrap intervals."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from mosaic import _paper_style
        _paper_style()
    except Exception:
        pass
    base, ci, n = shares_with_ci(path)
    order = base.sort_values().index
    vals = base[order] * 100
    lo = np.array([vals[t] - ci.get(t, (base[t], base[t]))[0] * 100 for t in order])
    hi = np.array([ci.get(t, (base[t], base[t]))[1] * 100 - vals[t] for t in order])
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#999999" if t == "Residual" else "#3B6EA8" for t in order]
    ax.barh([LABELS.get(t, t) for t in order], vals, color=colors,
            xerr=[np.clip(lo, 0, None), np.clip(hi, 0, None)], capsize=4)
    for i, v in enumerate(vals):
        ax.text(v + 1, i + 0.25, f"{v:.0f}%", fontsize=12)
    ax.set_xlabel("Share of variance in domain count (%)")
    ax.set_xlim(0, 100)
    ax.set_title(title or f"Categorical Type-II ANOVA ({n} pipelines)")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--figure":
        plot(args[1], args[2])
        sys.exit(0)
    for path in args:
        base, ci, n = shares_with_ci(path)
        print(path, n)
        for t in base.sort_values(ascending=False).index:
            lo, hi = ci.get(t, (np.nan, np.nan))
            print(f"  {t:10s} {base[t]*100:5.1f}%  [{lo*100:.1f}, {hi*100:.1f}]")
