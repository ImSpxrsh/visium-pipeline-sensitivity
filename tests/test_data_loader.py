"""
Unit tests for src/data_loader.py.

The load tests are marked ``network`` because they download the Visium demo.
Run only the offline tests with:  pytest -m "not network"
"""

import anndata as ad
import numpy as np
import pytest

from src.data_loader import prepare_data, validate_dataset


def _tiny_adata():
    x = np.random.RandomState(0).poisson(1.0, size=(5, 4)).astype(float)
    a = ad.AnnData(x)
    a.obs_names = [f"spot{i}" for i in range(5)]
    a.var_names = ["g0", "g1", "g1", "g2"]  # duplicate on purpose
    a.obsm["spatial"] = np.random.RandomState(1).rand(5, 2)
    return a


def test_validate_passes_on_well_formed_data():
    assert validate_dataset(_tiny_adata(), verbose=False) is True


def test_validate_fails_without_spatial():
    a = _tiny_adata()
    del a.obsm["spatial"]
    assert validate_dataset(a, verbose=False) is False


def test_prepare_makes_unique_names_and_stores_counts():
    prepared = prepare_data(_tiny_adata(), verbose=False)
    assert len(prepared.var_names) == len(set(prepared.var_names))
    assert "counts" in prepared.layers


@pytest.mark.network
def test_load_visium_hne():
    from src.data_loader import load_dataset

    adata, image = load_dataset("visium_hne", verbose=False)
    assert adata.n_obs > 0 and adata.n_vars > 0
    assert image is not None
    assert validate_dataset(adata, verbose=False) is True
