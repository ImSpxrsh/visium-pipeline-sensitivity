"""Utility functions: config loading, logging, file I/O, hashing, grid expansion."""

import hashlib
import json
import logging
from itertools import product
from pathlib import Path
from typing import Any, Dict, List

import yaml

logger = logging.getLogger(__name__)


def load_yaml(path: Path) -> Dict[str, Any]:
    """Load a YAML config file into a dict."""
    with open(path, "r") as f:
        return yaml.safe_load(f)


def save_json(data: Dict[str, Any], path: Path) -> None:
    """Save a dict as pretty-printed JSON, creating parent dirs as needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)


def load_json(path: Path) -> Dict[str, Any]:
    """Load a JSON file into a dict."""
    with open(path, "r") as f:
        return json.load(f)


def dict_to_hash(d: Dict[str, Any]) -> str:
    """
    Create a short deterministic hash from a dict (for naming spec directories).

    Example:
        dict_to_hash({"norm": "log1p", "resolution": 0.5}) -> "a7f3c2e1"
    """
    json_str = json.dumps(d, sort_keys=True, default=str)
    return hashlib.md5(json_str.encode()).hexdigest()[:8]


def setup_logging(log_file: Path, level: int = logging.INFO) -> None:
    """
    Configure the root logger with both file and console handlers.

    Idempotent: clears existing handlers so repeated calls don't duplicate lines.
    """
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Avoid duplicate handlers when called more than once in the same process.
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)


def expand_grid(spec: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    """
    Expand a dict of lists into all combinations (Cartesian product).

    Example:
        expand_grid({"a": [1, 2], "b": [10, 20]})
        -> [{"a": 1, "b": 10}, {"a": 1, "b": 20},
            {"a": 2, "b": 10}, {"a": 2, "b": 20}]

    Scalar (non-list) values are treated as a single-element axis.
    """
    keys = list(spec.keys())
    values = [spec[k] if isinstance(spec[k], list) else [spec[k]] for k in keys]
    return [dict(zip(keys, combo)) for combo in product(*values)]


def filter_specs_by_count(
    specs: List[Dict[str, Any]], max_count: int = 40
) -> List[Dict[str, Any]]:
    """Cap the number of specs to keep the sweep interpretable."""
    if len(specs) > max_count:
        logger.warning(
            "Too many specs (%d); keeping the first %d", len(specs), max_count
        )
        return specs[:max_count]
    return specs
