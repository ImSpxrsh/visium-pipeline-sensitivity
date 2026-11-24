"""Surrounding-tumor S3 gaps on the breast-cancer section, with spot-bootstrap CIs.

Writes results/breast_cancer/aggregate/pathology_margin_vs_invasive.json. The
bootstrap resamples spots independently, so it ignores spatial dependence and
the intervals are likely too narrow (as stated in the paper).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent
AGG = PROJECT / "results" / "breast_cancer" / "aggregate"
MARGIN, INVASIVE = "Surrounding tumor", "Invasive"


def main(n_boot: int = 1000, seed: int = 0) -> None:
    ps = pd.read_csv(AGG / "per_spot_stability.csv").set_index("spot_id")
    ann = pd.read_csv(PROJECT / "data" / "annotations" / "brca1_metadata.tsv",
                      sep="\t").set_index("ID")
    df = ps.join(ann, how="inner")
    s3, lab = df["s3_score"].values, df["annot_type"].values
    m, inv, rest = s3[lab == MARGIN], s3[lab == INVASIVE], s3[lab != MARGIN]
    rng = np.random.default_rng(seed)

    def boot(a, b):
        d = [rng.choice(b, len(b)).mean() - rng.choice(a, len(a)).mean() for _ in range(n_boot)]
        return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]

    out = {
        "margin_vs_invasive_gap": float(inv.mean() - m.mean()),
        "margin_vs_invasive_ci95": boot(m, inv),
        "margin_vs_rest_gap": float(rest.mean() - m.mean()),
        "margin_vs_rest_ci95": boot(m, rest),
        "n_margin": int(len(m)), "n_invasive": int(len(inv)), "n_rest": int(len(rest)),
        "n_boot": n_boot,
    }
    (AGG / "pathology_margin_vs_invasive.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
