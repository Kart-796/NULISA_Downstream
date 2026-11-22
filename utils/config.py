"""Shared configuration loader for Python scripts and notebooks.

    from utils.config import cfg, cfg_path, read_samples
    OUTPUT_BASE = cfg_path("multivelo_macs3_output")

Settings live in config/config.yaml (git-ignored; copy config/config.example.yaml).
Set PROJECT_CONFIG to use a different file. Relative paths are resolved against
the repository root, so scripts work no matter which folder they are started from.
"""
import csv
import os
from pathlib import Path

import yaml

REPO_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))
CONFIG_FILE = Path(os.environ.get("PROJECT_CONFIG", REPO_ROOT / "config" / "config.yaml"))

if not CONFIG_FILE.exists():
    raise FileNotFoundError(
        f"Config file not found: {CONFIG_FILE}\n"
        "Copy config/config.example.yaml to config/config.yaml and edit it."
    )

with open(CONFIG_FILE) as fh:
    cfg = yaml.safe_load(fh)


def _resolve(p):
    p = Path(os.path.expanduser(str(p)))
    return p if p.is_absolute() else (REPO_ROOT / p)


def cfg_path(key, create=False):
    """Absolute path (as str) for a key under `paths:` in the config."""
    value = (cfg.get("paths") or {}).get(key)
    if not value:
        raise KeyError(f"Missing 'paths: {key}' in {CONFIG_FILE}")
    p = _resolve(value)
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return str(p)


def read_samples(key):
    """Sample sheet for a key under `samples:` as a list of dicts."""
    f = (cfg.get("samples") or {}).get(key)
    if not f:
        raise KeyError(f"Missing 'samples: {key}' in {CONFIG_FILE}")
    with open(_resolve(f), newline="") as fh:
        return list(csv.DictReader(fh))


def sample_ids(key):
    return [row["sample_id"] for row in read_samples(key)]
