"""
Tests for the data-preparation helpers in src/visualization.py.

These exercise the SVG-concordance logic (no figure rendering needed). The
module uses a headless Agg backend, so importing it is safe in CI.
"""

import numpy as np
import pandas as pd

from src.visualization import _significant, svg_concordance_matrix


def _run(name, svg_df):
    return {"name": name, "svg_df": svg_df}


def test_significant_prefers_fdr_column():
    df = pd.DataFrame(
        {"I": [0.9, 0.1], "pval_norm": [0.9, 0.9], "pval_norm_fdr_bh": [0.01, 0.9]},
        index=["g0", "g1"],
    )
    sig = _significant(df)
    assert bool(sig["g0"]) is True
    assert bool(sig["g1"]) is False


def test_significant_falls_back_to_morans_i():
    df = pd.DataFrame({"I": [0.9, 0.1]}, index=["g0", "g1"])
    sig = _significant(df)
    assert bool(sig["g0"]) is True
    assert bool(sig["g1"]) is False


def _svg(genes, sig):
    return pd.DataFrame(
        {"I": np.linspace(1, 0, len(genes)),
         "pval_norm_fdr_bh": [0.01 if s else 0.9 for s in sig]},
        index=genes,
    )


def test_concordance_union_treats_absent_as_not_called():
    runs = [
        _run("s0", _svg(["g0", "g1"], [True, True])),
        _run("s1", _svg(["g0", "g2"], [True, True])),  # g1 absent, g2 present
    ]
    mat = svg_concordance_matrix(runs, common_only=False)
    assert set(mat.index) == {"g0", "g1", "g2"}
    # g1 untested in s1 -> counted as not called (False), not NaN.
    assert mat.loc["g1", "s1"] == False  # noqa: E712


def test_concordance_common_only_keeps_shared_genes():
    runs = [
        _run("s0", _svg(["g0", "g1"], [True, False])),
        _run("s1", _svg(["g0", "g2"], [False, True])),
    ]
    mat = svg_concordance_matrix(runs, common_only=True)
    assert list(mat.index) == ["g0"]  # only gene tested in both specs
    assert mat.loc["g0", "s0"] == True  # noqa: E712
    assert mat.loc["g0", "s1"] == False  # noqa: E712
