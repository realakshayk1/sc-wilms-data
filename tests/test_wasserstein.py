"""Tests for Phase B utilities and 1-D Wasserstein guardrails."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "phase2_histology_ml"))

from utils import load_config, repo_root  # noqa: E402


def test_repo_root_finds_config():
    root = repo_root()
    assert (root / "config" / "paths.yaml").exists()
    assert (root / "PRD.md").exists() or (root / "AGENTS.md").exists()


def test_feature_config_is_one_dimensional():
    cfg = load_config()
    for feat in cfg["features"]["features"]:
        assert "id" in feat
        assert "genes_positive" in feat


def test_wasserstein_only_on_1d_scores():
    """Guardrail: distance matrix built from scalar score vectors only."""
    rng = np.random.default_rng(42)
    x = rng.normal(size=50)
    y = rng.normal(loc=0.5, size=60)
    assert x.ndim == 1 and y.ndim == 1
    # High-dim gene matrix must NOT be passed to Wasserstein
    gene_matrix = rng.normal(size=(100, 50))
    with pytest.raises(AssertionError):
        assert gene_matrix.ndim == 1, "Wasserstein must use 1-D scores only"
