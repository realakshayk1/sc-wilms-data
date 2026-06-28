"""Shared utilities for Phase B (histology ML → ABM)."""

from __future__ import annotations

import json
import logging
import os
import random
from pathlib import Path
from typing import Any

import numpy as np
import yaml

LOG = logging.getLogger("wilms_abm")


def repo_root() -> Path:
    here = Path.cwd().resolve()
    for candidate in (here, here.parent, here.parent.parent):
        if (candidate / "config" / "paths.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate repo root (config/paths.yaml)")


def load_config() -> dict[str, Any]:
    root = repo_root()
    with open(root / "config" / "paths.yaml") as f:
        paths = yaml.safe_load(f)
    with open(root / "config" / "features.yaml") as f:
        features = yaml.safe_load(f)
    phase_b_path = root / "config" / "phase_b.yaml"
    phase_b = yaml.safe_load(phase_b_path.read_text()) if phase_b_path.exists() else {}
    return {"root": root, "paths": paths, "features": features, "phase_b": phase_b}


def resolve_path(cfg: dict[str, Any], rel: str) -> Path:
    return cfg["root"] / rel


def set_seed_logged(seed: int, label: str = "global") -> int:
    random.seed(seed)
    np.random.seed(seed)
    LOG.info("[seed] %s = %d", label, seed)
    return seed


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def is_demo_mode() -> bool:
    return "--demo" in os.sys.argv or os.environ.get("WILMS_DEMO", "0") == "1"


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )


def save_json(path: Path, obj: Any) -> None:
    ensure_dir(path.parent)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
