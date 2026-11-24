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


if __name__ == "__main__":
    for path in sys.argv[1:]:
        base, ci, n = shares_with_ci(path)
        print(path, n)
        for t in base.sort_values(ascending=False).index:
            lo, hi = ci.get(t, (np.nan, np.nan))
            print(f"  {t:10s} {base[t]*100:5.1f}%  [{lo*100:.1f}, {hi*100:.1f}]")
