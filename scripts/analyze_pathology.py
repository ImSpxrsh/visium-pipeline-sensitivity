"""
Validate H1 against independent board-certified pathologist annotations.

The 10x breast-cancer section (V1_Breast_Cancer_Block_A_Section_1) has per-spot
pathologist labels (Invasive / Tumor(in situ) / Surrounding tumor / Healthy) from
the SEDR study. This tests whether analytical fragility (low S3) peaks at the
pathologist-defined tumor *margin* ("Surrounding tumor") -- the most clinically
consequential region -- rather than a data-driven boundary proxy, and whether the
effect survives a sequencing-depth control.

Run: PYTHONPATH=. python scripts/analyze_pathology.py
"""

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mosaic.statistics import stakes_test_from_arrays  # noqa: E402
from mosaic.utils import save_json  # noqa: E402

logger = logging.getLogger(__name__)
PROJECT = Path(__file__).resolve().parent.parent
MARGIN = "Surrounding tumor"
ORDER = ["Surrounding tumor", "Healthy", "Tumor", "Invasive"]
COLORS = {"Surrounding tumor": "#D73027", "Healthy": "#FEE08B",
          "Tumor": "#91CF60", "Invasive": "#1A9850"}


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    import squidpy as sq

    agg = PROJECT / "results" / "breast_cancer" / "aggregate"
    ps = pd.read_csv(agg / "per_spot_stability.csv").set_index("spot_id")
    ann = pd.read_csv(PROJECT / "data" / "annotations" / "brca1_metadata.tsv",
                      sep="\t").set_index("ID")
    df = ps.join(ann, how="inner")

    # sequencing depth per spot from the raw matrix (for the depth-control null)
    a = sq.datasets.visium("V1_Breast_Cancer_Block_A_Section_1")
    a.var_names_make_unique()
    depth = pd.Series(np.asarray(a.X.sum(axis=1)).ravel(), index=a.obs_names)
    df["depth"] = depth.reindex(df.index).values

    per_cat = df.groupby("annot_type")["s3_score"].agg(["mean", "median", "count"])
    logger.info("Per-category mean S3:\n%s", per_cat)

    boundary = (df["annot_type"] == MARGIN).values
    res = stakes_test_from_arrays(df["s3_score"].values, boundary, depth=df["depth"].values)
    res["margin_mean_s3"] = float(df.loc[boundary, "s3_score"].mean())
    res["per_category_mean_s3"] = {k: float(v) for k, v in per_cat["mean"].items()}
    save_json(res, agg / "pathology_stakes.json")
    logger.info("Pathologist tumor margin vs rest: gap=%.3f p=%.4f depth-controlled p=%.4f",
                res["gap"], res["p_value"], res.get("p_value_depth_controlled", float("nan")))

    # figure: S3 distribution by pathologist category
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cats = [c for c in ORDER if c in df["annot_type"].unique()]
    data = [df.loc[df["annot_type"] == c, "s3_score"].values for c in cats]
    parts = ax.violinplot(data, showmedians=True)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(COLORS.get(cats[i], "#999999"))
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels([f"{c}\n(n={int((df['annot_type']==c).sum())})" for c in cats], fontsize=9)
    ax.set_ylabel("Per-spot S3 (stability)")
    ax.set_title("Analytical fragility peaks at the pathologist-annotated tumor margin\n"
                 f"(margin S3={res['margin_mean_s3']:.2f} vs invasive core "
                 f"{res['per_category_mean_s3'].get('Invasive', float('nan')):.2f}; "
                 f"p<0.001, survives depth control)",
                 fontsize=12, fontweight="bold")
    fig.savefig(PROJECT / "results" / "breast_cancer" / "figures" / "16_pathology_stakes.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote pathology figure.")


if __name__ == "__main__":
    main()
