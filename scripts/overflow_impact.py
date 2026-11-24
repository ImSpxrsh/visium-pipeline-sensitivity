"""Compare per-spot S3 before and after the overflow-label fix on a finished sweep.

Usage: python scripts/overflow_impact.py results/full_hne [results/full_breast_cancer ...]
Reads only obs['leiden'] from each spec, so it is fast and light on memory.
"""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mosaic.stability import _align_labels_to_reference  # noqa: E402
from mosaic.statistics import _boundary_flags  # noqa: E402


def old_align(labels, ref_labels):
    out = _align_labels_to_reference(labels, ref_labels, spec_index=0)
    n_ref = int(pd.Categorical(ref_labels).codes.max()) + 1
    out[out >= n_ref] = n_ref  # the old shared overflow label
    return out


def s3(mat):
    vals = []
    for row in mat:
        _, c = np.unique(row, return_counts=True)
        vals.append(c.max() / len(row))
    return np.array(vals)


for root in sys.argv[1:]:
    runs = sorted(Path(root, "runs").glob("spec_*/adata_clustered.h5ad"),
                  key=lambda p: int(p.parent.name.split("_")[1]))
    labs, coords = [], None
    for p in runs:
        a = ad.read_h5ad(p, backed="r")
        labs.append(np.asarray(a.obs["leiden"].values))
        if coords is None:
            coords = np.asarray(a.obsm["spatial"])
        a.file.close()
    ref = labs[0]
    old = np.array([old_align(l, ref) for l in labs]).T
    new = np.array([_align_labels_to_reference(l, ref, spec_index=i)
                    for i, l in enumerate(labs)]).T
    n_ref = int(pd.Categorical(ref).codes.max()) + 1
    s_old, s_new = s3(old), s3(new)
    out = {"root": root, "n_specs": len(labs),
           "overflow_entry_frac": float((old == n_ref).mean()),
           "spots_changed_frac": float((np.abs(s_old - s_new) > 1e-12).mean())}
    for tag, s, mat in (("old", s_old, old), ("new", s_new, new)):
        modal = np.array([np.unique(r, return_counts=True)[0][np.unique(r, return_counts=True)[1].argmax()] for r in mat])
        b = _boundary_flags(coords, modal, k=6)
        out[tag] = {"mean": float(s.mean()), "frac_ge_0.8": float((s >= 0.8).mean()),
                    "frac_lt_0.5": float((s < 0.5).mean()),
                    "boundary": float(s[b].mean()), "interior": float(s[~b].mean()),
                    "gap": float(s[~b].mean() - s[b].mean())}
    print(json.dumps(out, indent=1))
