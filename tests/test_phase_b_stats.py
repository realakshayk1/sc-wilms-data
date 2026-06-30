"""Tests for Phase B validation statistics (phase_b_stats.py): DeLong CI,
permutation AUC p-value, and paired DeLong."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
PHASE_B = ROOT / "phase2_histology_ml"
sys.path.insert(0, str(PHASE_B))


def _load():
    path = PHASE_B / "phase_b_stats.py"
    if not path.exists():
        pytest.skip("phase_b_stats.py missing")
    spec = importlib.util.spec_from_file_location("phase_b_stats", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:  # missing scipy/sklearn
        pytest.skip(f"cannot import: {e}")
    return mod


def test_delong_auc_matches_sklearn():
    """DeLong AUC point estimate must equal the Mann-Whitney/sklearn AUC."""
    sk = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 200)
    score = rng.normal(y * 0.8, 1.0)  # informative but noisy
    out = _load().delong_auc_ci(y, score)
    assert abs(out["auc"] - sk.roc_auc_score(y, score)) < 1e-9
    # CI brackets the point estimate and stays in [0, 1]
    assert 0.0 <= out["ci_low"] <= out["auc"] <= out["ci_high"] <= 1.0


def test_delong_perfect_and_chance():
    m = _load()
    # perfect separation -> AUC 1.0
    y = np.array([0, 0, 0, 1, 1, 1])
    perfect = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
    assert m.delong_auc_ci(y, perfect)["auc"] == 1.0
    # a wider, noisier sample's CI should be tighter than a tiny one's
    rng = np.random.default_rng(1)
    big_y = rng.integers(0, 2, 400); big_s = rng.normal(big_y, 1)
    small_y = big_y[:20]; small_s = big_s[:20]
    assert m.delong_auc_ci(big_y, big_s)["se"] < m.delong_auc_ci(small_y, small_s)["se"]


def test_permutation_p_separates_signal_from_noise():
    m = _load()
    y = np.r_[np.zeros(25), np.ones(25)].astype(int)
    strong = np.r_[np.zeros(25), np.ones(25)] + 0.01  # near-perfect ranking
    p_strong = m.permutation_auc_p(y, strong, n_perm=2000, seed=0)["p_value"]
    assert p_strong < 0.01
    rng = np.random.default_rng(2)
    noise = rng.normal(size=50)  # unrelated to y
    p_noise = m.permutation_auc_p(y, noise, n_perm=2000, seed=0)["p_value"]
    assert p_noise > 0.05


def test_paired_delong_direction_and_symmetry():
    m = _load()
    rng = np.random.default_rng(3)
    y = rng.integers(0, 2, 300)
    good = rng.normal(y * 1.2, 1.0)   # stronger classifier
    weak = rng.normal(y * 0.2, 1.0)   # weaker classifier
    res = m.delong_paired_test(y, good, weak)
    assert res["auc_a"] > res["auc_b"]          # A is the better model
    assert res["delta"] > 0
    # identical inputs -> zero difference, non-significant
    same = m.delong_paired_test(y, good, good.copy())
    assert abs(same["delta"]) < 1e-9
    assert same["p_value"] > 0.99 or np.isnan(same["p_value"])
