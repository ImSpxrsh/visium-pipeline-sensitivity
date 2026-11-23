"""Offline test for the single-spec option used by cluster job arrays."""

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np


def _load_main():
    path = Path(__file__).resolve().parent.parent / "scripts" / "main.py"
    spec = importlib.util.spec_from_file_location("mosaic_main", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_only_index_runs_one_spec(tmp_path, monkeypatch):
    main = _load_main()
    calls = []

    def fake_pipeline(adata, spec, random_state=42, verbose=False):
        calls.append(spec["name"])
        return {"adata": adata, "svg_results": None, "metadata": {"spec": spec}}

    monkeypatch.setattr(main, "run_pipeline", fake_pipeline)
    adata = ad.AnnData(np.ones((5, 3), dtype=np.float32))
    specs = [{"name": f"s{i}", "x": i} for i in range(4)]

    out = main.run_multiverse_sweep(adata, specs, tmp_path, verbose=False, only_index=2)

    assert calls == ["s2"]
    assert len(out) == 1
    dirs = [p.name for p in tmp_path.glob("spec_*")]
    assert len(dirs) == 1 and dirs[0].startswith("spec_02_")
    assert (tmp_path / "failure_manifest_spec2.json").exists()

    # Running all specs afterwards reuses the cached one and runs the other three.
    calls.clear()
    out = main.run_multiverse_sweep(adata, specs, tmp_path, verbose=False)
    assert sorted(calls) == ["s0", "s1", "s3"]
    assert len(out) == 4


def test_lean_aggregation_matches_full(tmp_path):
    """Slimmed AnnData gives the same stability outputs as the full object."""
    import pandas as pd

    from mosaic.stability import StabilityScorer

    main = _load_main()
    rng = np.random.default_rng(0)
    n = 40
    full = []
    for s in range(5):
        a = ad.AnnData(rng.random((n, 10)).astype(np.float32))
        a.obs["leiden"] = pd.Categorical(rng.integers(0, 4, n).astype(str))
        a.obs["cluster"] = pd.Categorical(rng.integers(0, 3, n).astype(str))
        a.obsm["spatial"] = rng.random((n, 2))
        a.layers["counts"] = a.X.copy()
        a.uns["spatial"] = {"img": np.zeros((50, 50, 3))}
        a.uns["svg_results"] = pd.DataFrame(
            {"I": rng.random(10), "pval_norm_fdr_bh": rng.random(10)},
            index=[f"g{i}" for i in range(10)],
        )
        full.append(a)

    def outputs(adatas):
        sc = StabilityScorer(verbose=False)
        for i, a in enumerate(adatas):
            sc.add_result({"name": str(i)}, a, {"n_clusters_found": i})
        return (
            sc.compute_per_spot_stability("leiden"),
            sc.compute_per_feature_stability(),
            sc.compute_partition_concordance("leiden"),
            sc.concordance_with_reference("leiden", "cluster"),
        )

    f = outputs(full)
    s = outputs([main._slim(a) for a in full])
    pd.testing.assert_frame_equal(f[0], s[0])
    pd.testing.assert_frame_equal(f[1], s[1])
    assert f[2] == s[2] and f[3] == s[3]
    slim = main._slim(full[0])
    assert slim.X is None
    assert list(slim.uns.keys()) == ["svg_results"]
