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


def test_rules_half_max_per_tumor_and_fallback():
    """The oxygen-necrosis and pressure half-maxes come from the per-tumor `half_max` block
    (17_positives_to_abm.py); absent -> documented defaults."""
    mod = _load("03_emit_rules.py")
    rates = {"proliferation_rate": 0.05, "apoptosis_rate": 0.001, "adhesion_strength": 0.3}
    # fallback (no half_max) uses DEFAULT_HM
    base = {r["behavior"] + "|" + r["signal"]: r["half_max"] for r in mod.rules_for("blastemal", rates)}
    assert base["necrosis|oxygen"] == pytest.approx(mod.DEFAULT_HM["oxygen_necrosis"])
    assert base["cycle entry|pressure"] == pytest.approx(mod.DEFAULT_HM["pressure_cycle"])
    # per-tumor override flows into the emitted rows
    hm = {"oxygen_cycle": 12.0, "oxygen_necrosis": 3.5, "pressure_cycle": 0.7}
    got = {r["behavior"] + "|" + r["signal"]: r["half_max"] for r in mod.rules_for("blastemal", rates, hm)}
    assert got["necrosis|oxygen"] == pytest.approx(3.5)
    assert got["cycle entry|pressure"] == pytest.approx(0.7)
    assert got["cycle entry|oxygen"] == pytest.approx(12.0)


def test_rules_igf2_ecm_gated_by_substrates():
    """v1.1 IGF2 / ECM rules appear only when those substrates are configured, and layer
    additively (base 3 rules -> 5 with both)."""
    mod = _load("03_emit_rules.py")
    rates = {"proliferation_rate": 0.05, "apoptosis_rate": 0.001, "adhesion_strength": 0.3}
    base = mod.rules_for("blastemal", rates)
    assert len(base) == 3                                  # no substrate set -> oxygen/pressure only
    full = mod.rules_for("blastemal", rates, None, {"IGF2", "ECM"})
    sigs = {(r["signal"], r["behavior"]) for r in full}
    assert ("IGF2", "cycle entry") in sigs
    assert ("ECM", "migration speed") in sigs
    igf = [r for r in full if r["signal"] == "IGF2"][0]
    assert igf["response"] == "increases"
    assert igf["saturation_value"] == pytest.approx(1.5 * 0.05)   # anchored to base proliferation
    ecm = [r for r in full if r["signal"] == "ECM"][0]
    assert ecm["response"] == "decreases" and float(ecm["saturation_value"]) == 0.0


def test_build_xml_substrates_and_secretion():
    """04 writes IGF2/ECM microenvironment variables and a per-cell secretion block with an
    entry per substrate (oxygen + configured extras)."""
    mod = _load("04_build_model.py")
    tumor = {"high_grade_regime": False,
             "cell_types": {c: {"proliferation_rate": 0.05, "apoptosis_rate": 0.001,
                                "adhesion_strength": 0.3, "migration_speed": 0.4,
                                "igf_uptake_rate": 0.001, "ecm_secretion_rate":
                                (0.001 if c == "stromal" else 0.0)}
                            for c in mod.COMPARTMENTS}}
    substrates = {"IGF2": {"diffusion": 1000.0, "decay": 0.01, "initial": 1.0,
                           "boundary": 1.0, "dirichlet": True},
                  "ECM": {"diffusion": 0.0, "decay": 0.0, "initial": 0.0,
                          "boundary": 0.0, "dirichlet": False}}
    root = mod.build_xml("S1", tumor, {"x_max": 500, "y_max": 500},
                         {"max_time_min": 100, "save_interval_min": 10}, substrates)
    names = {v.get("name") for v in root.iter("variable")}
    assert {"oxygen", "IGF2", "ECM"}.issubset(names)
    # each cell definition has a secretion entry per substrate (oxygen + IGF2 + ECM = 3)
    for cd in root.iter("cell_definition"):
        subs = [s.get("name") for s in cd.iter("substrate")]
        assert subs == ["oxygen", "IGF2", "ECM"]
    # stromal secretes ECM; tumor cells do not
    for cd in root.iter("cell_definition"):
        for s in cd.iter("substrate"):
            if s.get("name") == "ECM":
                rate = float(s.find("secretion_rate").text)
                assert (rate > 0) == (cd.get("name") == "stromal")


def test_clustering_index_segregation_vs_mixing():
    v = _load("07_validate.py")
    cats = ["blastemal", "stromal"]
    # two separated homotypic blocks -> clustering index > 1 for each type
    left = np.column_stack([np.zeros(40), np.linspace(0, 40, 40)])
    right = np.column_stack([np.full(40, 100.0), np.linspace(0, 40, 40)])
    seg = v.clustering_index(np.vstack([left, right]),
                             np.array(["blastemal"] * 40 + ["stromal"] * 40), cats, k=4)
    assert (seg["clustering_index"] > 1.2).all()
    # interleaved checkerboard -> index near 1 (well mixed)
    xs, ys = np.meshgrid(np.arange(10), np.arange(10))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    lab = np.where((coords[:, 0].astype(int) + coords[:, 1].astype(int)) % 2 == 0,
                   "blastemal", "stromal")
    mix = v.clustering_index(coords, lab, cats, k=4)
    assert (mix["clustering_index"] < 1.05).all()


def test_radial_invasiveness_reacts_to_outliers():
    """Against a FIXED reference radius (as the sim uses the t0 median at later timepoints),
    a compact mass has no projections beyond it; adding far spokes creates invasive ones."""
    v = _load("07_validate.py")
    rng = np.random.default_rng(0)
    disc = rng.normal(0, 3, (200, 2))                     # p95 radius ~ 7 um
    ref = 15.0                                            # fixed reference the disc never reaches
    compact = v.radial_invasiveness(disc, reference_radius=ref)
    assert compact["n_invasive_projections"] == 0
    assert compact["invasive_fraction"] == 0.0
    spokes = np.array([[80.0, 0.0], [-80.0, 0.0], [0.0, 80.0], [0.0, -80.0]])
    invasive = v.radial_invasiveness(np.vstack([disc, spokes]), reference_radius=ref)
    assert invasive["n_invasive_projections"] >= 4        # the four spokes' sectors
    assert invasive["invasive_fraction"] > 0.0
    assert invasive["radial_p95_over_ref"] > compact["radial_p95_over_ref"]


def test_zscore_col_optional_and_neutral():
    import pandas as pd
    import importlib.util as _il
    path = ROOT / "phase2_histology_ml" / "17_positives_to_abm.py"
    sys.path.insert(0, str(ROOT / "phase2_histology_ml"))
    spec = _il.spec_from_file_location("positives_to_abm", path)
    try:
        mod = _il.module_from_spec(spec)
        spec.loader.exec_module(mod)
    except Exception as e:  # pragma: no cover
        pytest.skip(f"cannot import 17_positives_to_abm: {e}")
    df = pd.DataFrame({"present": [1.0, 2.0, 3.0, 4.0]})
    # absent column -> all zeros (neutral)
    assert (mod.zscore_col(df, "missing") == 0.0).all()
    # present column -> proper z-score (mean 0, unit population std)
    z = mod.zscore_col(df, "present")
    assert z.mean() == pytest.approx(0.0, abs=1e-9)
    assert z.std(ddof=0) == pytest.approx(1.0, abs=1e-9)


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


def test_cooccurrence_detects_colocation_and_segregation():
    v = _load("07_validate.py")
    # two spatially separated blocks -> at short range same-type co-locates (>1),
    # cross-type segregates (<1)
    left = np.column_stack([np.zeros(60), np.linspace(0, 60, 60)])
    right = np.column_stack([np.full(60, 300.0), np.linspace(0, 60, 60)])
    coords = np.vstack([left, right])
    labels = np.array(["blastemal"] * 60 + ["stromal"] * 60)
    df = v.co_occurrence(coords, labels, ["blastemal", "stromal"],
                         radii=[10, 30, 60, 120, 240, 480])
    near = df[df.r_um <= 60]
    same = near[(near.cond == "blastemal") & (near.exp == "blastemal")]["cooccur"].mean()
    cross = near[(near.cond == "blastemal") & (near.exp == "stromal")]["cooccur"].mean()
    assert same > 1.0 > cross


def test_squidpy_extras_if_installed():
    """When squidpy is present, the optional backend returns a coherent nhood z-matrix
    and a Ripley's L table. Skipped otherwise (squidpy is an optional dependency)."""
    v = _load("07_validate.py")
    if not v.has_squidpy():
        pytest.skip("squidpy not installed (optional)")
    rng = np.random.default_rng(0)
    left = np.column_stack([rng.normal(0, 5, 120), rng.normal(0, 5, 120)])
    right = np.column_stack([rng.normal(200, 5, 120), rng.normal(0, 5, 120)])
    coords = np.vstack([left, right])
    labels = np.array(["blastemal"] * 120 + ["stromal"] * 120)
    ex = v.squidpy_extras(coords, labels, v.COMPARTMENTS, knn=6, n_perm=100, seed=0)
    z = ex["nhood_z"]
    assert "blastemal" in z.index and "stromal" in z.index
    # same-type adjacency enriched vs cross-type (segregated blocks)
    assert z.loc["blastemal", "blastemal"] > z.loc["blastemal", "stromal"]
    assert not ex["ripley_L"].empty


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
