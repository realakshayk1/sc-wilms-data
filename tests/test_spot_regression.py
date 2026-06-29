"""Tests for Phase B spot-level composition regression helpers."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
PHASE_B = ROOT / "phase2_histology_ml"
sys.path.insert(0, str(PHASE_B))


def _load_module():
    path = PHASE_B / "12_spot_composition_regression.py"
    if not path.exists():
        pytest.skip("regression script missing")
    spec = importlib.util.spec_from_file_location("spot_reg", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:  # missing optional deps
        pytest.skip(f"cannot import: {e}")
    return mod


def test_zscore_softmax_target_is_a_valid_composition():
    mod = _load_module()
    rng = np.random.default_rng(0)
    n = 500
    sig = pd.DataFrame({
        "spot_id": [f"s{i}" for i in range(n)],
        "blastemal_program": rng.normal(2.0, 1, n),   # deliberately offset high
        "epithelial_program": rng.normal(-0.5, 1, n),
        "stromal_program": rng.normal(0.0, 1, n),
    })
    out = mod.zscore_softmax_target(sig)
    cols = [f"y_{s}" for s in mod.CELL_STATES]
    # rows are a probability simplex
    assert np.allclose(out[cols].sum(axis=1), 1.0, atol=1e-6)
    assert (out[cols].to_numpy() >= 0).all()
    # z-scoring removes the additive offset: no single compartment dominates every spot
    dom = out[cols].to_numpy().argmax(1)
    assert len(np.unique(dom)) == 3, "z-scored target should not collapse to one class"


def test_spot_morphology_aggregates_per_spot():
    mod = _load_module()
    nuc = pd.DataFrame({
        "spot_id": ["a", "a", "b"],
        "sample_id": ["S1", "S1", "S1"],
        "subdiagnosis": ["favorable"] * 3,
        **{m: [1.0, 3.0, 2.0] for m in mod.MORPH},
    })
    agg = mod.spot_morphology(nuc)
    assert set(agg["spot_id"]) == {"a", "b"}
    row_a = agg[agg["spot_id"] == "a"].iloc[0]
    assert row_a["n_nuclei"] == 2
    assert abs(row_a["area_mean"] - 2.0) < 1e-9
