import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "network: test downloads data and needs a network connection"
    )
