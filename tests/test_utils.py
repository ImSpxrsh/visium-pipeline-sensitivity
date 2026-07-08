"""Unit tests for src/utils.py (no network required)."""

from src.utils import dict_to_hash, expand_grid, filter_specs_by_count


def test_dict_to_hash_is_deterministic_and_order_independent():
    a = {"norm": "log1p", "res": 0.5}
    b = {"res": 0.5, "norm": "log1p"}
    assert dict_to_hash(a) == dict_to_hash(b)
    assert len(dict_to_hash(a)) == 8


def test_dict_to_hash_differs_on_change():
    assert dict_to_hash({"res": 0.5}) != dict_to_hash({"res": 0.8})


def test_expand_grid_cartesian_product():
    grid = expand_grid({"a": [1, 2], "b": [10, 20]})
    assert len(grid) == 4
    assert {"a": 1, "b": 10} in grid
    assert {"a": 2, "b": 20} in grid


def test_expand_grid_scalar_axis():
    grid = expand_grid({"a": [1, 2], "b": 99})
    assert len(grid) == 2
    assert all(g["b"] == 99 for g in grid)


def test_filter_specs_caps_count():
    specs = [{"i": i} for i in range(50)]
    assert len(filter_specs_by_count(specs, max_count=30)) == 30
    assert len(filter_specs_by_count(specs[:10], max_count=30)) == 10
