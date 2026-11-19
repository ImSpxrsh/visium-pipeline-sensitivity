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
