"""Tests for Phase C ABM authoring (abm_utils + rules emitter + placement).

Covers the pure, data-independent invariants: the pixel->micron affine, the Visium
coordinate loader on a synthetic library, and the grammar-rule schema. A placement smoke
test runs only if the real processed data is present (skipped otherwise)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
PHASE_C = ROOT / "phase3_abm"
sys.path.insert(0, str(PHASE_C))


def _load(fname: str):
    path = PHASE_C / fname
    spec = importlib.util.spec_from_file_location(path.stem.lstrip("0123456789_") or "mod", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:  # pragma: no cover
        pytest.skip(f"cannot import {fname}: {e}")
    return mod


def test_um_per_pixel():
    au = _load("abm_utils.py")
    # 55 um spot spanning 110 px -> 0.5 um/px
    assert au.um_per_pixel({"spot_diameter_fullres": 110.0}, 55.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        au.um_per_pixel({"spot_diameter_fullres": 0.0}, 55.0)


def test_load_spot_coords_affine(tmp_path):
    au = _load("abm_utils.py")
    spatial = tmp_path / "spatial"
    spatial.mkdir(parents=True)
    # two spots 100 px apart in x; 0.5 um/px -> 50 um apart
    pd.DataFrame([
        ["AAA-1", 1, 0, 0, 100, 200],
        ["BBB-1", 0, 0, 1, 200, 200],
    ]).to_csv(spatial / "tissue_positions_list.csv", header=False, index=False)
    (spatial / "scalefactors_json.json").write_text(json.dumps({"spot_diameter_fullres": 110.0}))

    coords = au.load_spot_coords_um(tmp_path, spot_diameter_um=55.0)
    assert coords.attrs["um_per_pixel"] == pytest.approx(0.5)
    assert coords.loc["AAA-1", "x_um"] == pytest.approx(50.0)   # 100 px * 0.5
    assert coords.loc["BBB-1", "x_um"] == pytest.approx(100.0)  # 200 px * 0.5
    assert coords.loc["AAA-1", "y_um"] == pytest.approx(100.0)  # 200 px * 0.5
    assert int(coords.loc["AAA-1", "in_tissue"]) == 1


def test_rules_schema():
    mod = _load("03_emit_rules.py")
    rows = mod.rules_for("blastemal", {"proliferation_rate": 0.05,
                                       "apoptosis_rate": 0.001, "adhesion_strength": 0.3})
    assert len(rows) == 3
    for r in rows:
        for col in mod.HEADER:
            assert col in r                      # every PhysiCell column present
        assert r["response"] in ("increases", "decreases")
        assert float(r["saturation_value"]) >= 0.0
        assert int(r["apply_to_dead"]) in (0, 1)
        assert r["cell_type"] == "blastemal"
    # oxygen->cycle-entry saturation is 1.5x the base proliferation rate
    ox = [r for r in rows if r["signal"] == "oxygen" and r["behavior"] == "cycle entry"][0]
    assert ox["saturation_value"] == pytest.approx(1.5 * 0.05)


def test_placement_smoke_if_data_present():
    """Deterministic, in-bounds placement with fractions summing to 1 — real data only."""
    au = _load("abm_utils.py")
    cfg = au.load_config()
    sig_path = au.resolve_path(cfg, "data/processed/spot_signatures.parquet")
    spatial_root = au.resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"])
    if not sig_path.exists() or not spatial_root.exists():
        pytest.skip("processed Visium data not present")
    place = _load("02_place_agents.py")
    libs = au.discover_library_dirs(spatial_root)
    if not libs:
        pytest.skip("no Visium libraries on disk")
    sid = next(iter(libs))
    sig = pd.read_parquet(sig_path, columns=[
        "spot_id", "sample_id", "barcode", "in_tissue",
        *cfg["phase_c"]["deconvolution"]["frac_cols"]])
    seed = int(cfg["phase_c"]["seed"])
    a = place.place_tumor(cfg, sid, libs[sid], sig, au.rng_for(seed, "t"))
    b = place.place_tumor(cfg, sid, libs[sid], sig, au.rng_for(seed, "t"))
    if a is None:
        pytest.skip("no placeable spots for first library")
    assert list(a.columns) == ["x", "y", "z", "cell_type"]
    assert (a["x"] >= 0).all() and (a["y"] >= 0).all()
    assert set(a["cell_type"]).issubset(set(au.COMPARTMENTS))
    assert a.equals(b)                                    # deterministic given seed


def test_neighbor_enrichment_detects_segregation():
    v = _load("07_validate.py")
    # two spatially separated blocks: same-type neighbours enriched, cross-type depleted
    left = np.column_stack([np.zeros(50), np.linspace(0, 50, 50)])
    right = np.column_stack([np.full(50, 100.0), np.linspace(0, 50, 50)])
    coords = np.vstack([left, right])
    labels = np.array(["blastemal"] * 50 + ["stromal"] * 50)
    df = v.neighbor_enrichment(coords, labels, ["blastemal", "stromal"], k=4,
                               n_perm=300, seed=0)
    cross = df[(df.cat_a == "blastemal") & (df.cat_b == "stromal")]["z"].iloc[0]
    same = df[(df.cat_a == "blastemal") & (df.cat_b == "blastemal")]["z"].iloc[0]
    assert cross < 0 < same                               # segregation


def test_emergent_test_direction_and_fdr():
    import pandas as pd
    v = _load("07_validate.py")
    rng = np.random.default_rng(0)
    # 20 tumors; QoI 'growth' higher in group 1, 'noise' unrelated
    grp = np.r_[np.ones(10), np.zeros(10)].astype(int)
    df = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(20)],
        "anaplastic": grp,
        "growth": np.r_[rng.normal(5, 1, 10), rng.normal(1, 1, 10)],
        "noise": rng.normal(0, 1, 20)})
    res = v.patient_level_emergent_test(df, "anaplastic", ["growth", "noise"]).set_index("qoi")
    assert res.loc["growth", "p_value"] < 0.05
    assert res.loc["growth", "cliffs_delta"] > 0.5        # group 1 larger
    assert res.loc["noise", "p_value"] > 0.05
    assert (res["p_bh"] >= res["p_value"]).all()          # BH >= raw p


def test_cliffs_delta_bounds():
    v = _load("07_validate.py")
    assert v.cliffs_delta([3, 4, 5], [0, 1, 2]) == 1.0
    assert v.cliffs_delta([0, 1, 2], [3, 4, 5]) == -1.0


def test_generate_sweep_design():
    mod = _load("05_uq.py")
    df = mod.generate_sweep("S1", {"a": 10.0, "b": 0.4}, [10, 20])
    # 1 baseline + 2 params * 2 levels * 2 signs
    assert len(df) == 1 + 2 * 2 * 2
    assert (df["run_id"] == "S1__base").sum() == 1
    a_p10 = df[(df.param == "a") & (df.pct == 10)].iloc[0]
    assert a_p10["perturbed_value"] == pytest.approx(11.0)
    a_m20 = df[(df.param == "a") & (df.pct == -20)].iloc[0]
    assert a_m20["perturbed_value"] == pytest.approx(8.0)


def test_compute_qoi():
    import pandas as pd
    q = _load("qoi_extract.py")
    census = pd.DataFrame({"time": [0, 60, 120], "count": [100, 150, 200],
                           "blastemal_count": [60, 90, 120]})
    pos = pd.DataFrame({"x": [0.0, 10.0, -10.0, 0.0], "y": [0.0, 0.0, 0.0, 20.0]})
    out = q.compute_qoi(census, pos)
    assert out["final_total_cells"] == 200
    assert out["fold_growth"] == pytest.approx(2.0)
    assert out["growth_auc"] == pytest.approx(150.0)      # trapezoid mean
    assert out["radial_extent_um"] > out["radius_median_um"]
    assert out["final_blastemal_frac"] == pytest.approx(0.6)
