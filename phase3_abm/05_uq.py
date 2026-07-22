#!/usr/bin/env python3
"""Stage 5: uncertainty-quantification / sensitivity sweep manifest (per tumor).

Emits a one-at-a-time (OAT) sensitivity design: each uncertain parameter is perturbed by
+/- each pct level while the others hold at base, plus one baseline run. The manifest lists
every run the cluster should execute; the run wrapper (06_run_cohort) reads it and applies
the override to that tumor's rules.csv / PhysiCell_settings.xml before launching.

The design generation is pure and unit-tested here; only the execution is cluster-side.

Writes results/abm/<sample_id>/uq/sweep_manifest.csv.

Usage: python 05_uq.py [--sample SCPCS000168 ...]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import ensure_dir, load_config, resolve_path, setup_logging  # noqa: E402

# bounded per-tumor transforms (mirrored from phase2_histology_ml/17_positives_to_abm.py so the
# virtual-cohort draws use the SAME intrinsic base-rate mapping as real-tumor seeding)
PROLIF_K, PROLIF_LO, PROLIF_HI = 0.60, 0.40, 2.50
APOP_K, APOP_LO, APOP_HI = 0.50, 0.40, 2.00


def _norm_cdf(z: np.ndarray) -> np.ndarray:
    """Standard-normal CDF without scipy (vectorized math.erf)."""
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))


def _empirical_quantile(col_sorted: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Linear-interpolated empirical inverse-CDF; preserves the lever's own (maybe bimodal) marginal."""
    grid = (np.arange(len(col_sorted)) + 0.5) / len(col_sorted)
    return np.interp(u, grid, col_sorted)


def generate_virtual_cohort(scores: pd.DataFrame, corr: pd.DataFrame, levers: list[str],
                            n_draws: int, seed: int, ridge: float = 0.10) -> pd.DataFrame:
    """Empirical-marginal Gaussian copula: draw n_draws correlated lever vectors that preserve
    (a) each lever's observed marginal (bimodality included) and (b) the shrunk coupling structure.

    Pure/deterministic given seed — unit-tested. Returns one row per synthetic tumor with the drawn
    lever z-scores (extrinsic axes + PhysiCell params are added by apply_transfer)."""
    R = corr.loc[levers, levers].to_numpy(dtype=float)
    R = (1.0 - ridge) * R + ridge * np.eye(len(levers))     # ridge-shrink -> positive definite
    L = np.linalg.cholesky(R)
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((n_draws, len(levers))) @ L.T
    U = _norm_cdf(Z)
    out = {}
    for j, lev in enumerate(levers):
        col = np.sort(scores[lev].dropna().to_numpy())
        out[lev] = _empirical_quantile(col, U[:, j])
    df = pd.DataFrame(out)
    df.insert(0, "draw_id", [f"vc{i:04d}" for i in range(n_draws)])
    return df


def apply_transfer(draws: pd.DataFrame, transfer: dict, o2p: dict,
                   base_hm: dict) -> pd.DataFrame:
    """Derive extrinsic axes from levers via the measured transfer, then map to PhysiCell params
    with the existing bounded transforms (config/phase_c.yaml:omics_to_params). Axes are FUNCTIONS
    of the coupled levers, so the couplings propagate into the sampled parameters."""
    d = draws.copy()

    def axis_z(axis: str) -> np.ndarray:
        betas = transfer.get("measured", {}).get(axis, {})
        z = np.zeros(len(d))
        for lever, spec in betas.items():
            if lever in d:
                z = z + float(spec["beta"]) * d[lever].to_numpy()
        return z

    crowd_z = axis_z("crowding_sensitivity")
    hypox_z = axis_z("hypoxia_tolerance")
    press_k = float(o2p.get("pressure_halfmax_k", 0.30))
    necr_k = float(o2p.get("necrosis_halfmax_k", 0.30))
    hm_lo, hm_hi = o2p.get("halfmax_bounds", [0.5, 1.5])
    emt_adh_k = float(o2p.get("emt_adhesion_k", 0.25))
    emt_mot_k = float(o2p.get("emt_motility_k", 0.40))
    adh_lo, adh_hi = o2p.get("adhesion_bounds", [0.4, 1.6])
    mot_lo, mot_hi = o2p.get("motility_bounds", [0.4, 2.0])
    igf_k = float(o2p.get("igf_uptake_k", 0.40))
    up_lo, up_hi = o2p.get("uptake_bounds", [0.4, 2.5])

    emt = d["emt_axis"].to_numpy() if "emt_axis" in d else np.zeros(len(d))
    d["crowding_z"] = np.round(crowd_z, 4)
    d["hypoxia_z"] = np.round(hypox_z, 4)
    d["proliferation_mult"] = np.clip(1 + PROLIF_K * d["proliferation"], PROLIF_LO, PROLIF_HI).round(4)
    d["apoptosis_mult"] = np.clip(1 + APOP_K * d["tp53_target"], APOP_LO, APOP_HI).round(4)
    d["adhesion_mult"] = np.clip(1 - emt_adh_k * emt, adh_lo, adh_hi).round(4)
    d["motility_mult"] = np.clip(1 + emt_mot_k * emt, mot_lo, mot_hi).round(4)
    d["igf_uptake_mult"] = np.clip(1 + igf_k * d["igf"], up_lo, up_hi).round(4)
    d["pressure_half_max"] = (float(base_hm["pressure_half_max"])
                              * np.clip(1 - press_k * crowd_z, hm_lo, hm_hi)).round(4)
    d["necrosis_half_max"] = (float(base_hm["oxygen_necrosis_half_max"])
                              * np.clip(1 - necr_k * hypox_z, hm_lo, hm_hi)).round(4)
    return d


def generate_sweep(sample_id: str, base_values: dict[str, float],
                   pct_levels: list[float]) -> pd.DataFrame:
    """OAT design: baseline + (param x +/-pct) rows. Deterministic, order-stable."""
    rows = [{"run_id": f"{sample_id}__base", "sample_id": sample_id, "param": "(baseline)",
             "pct": 0.0, "factor": 1.0, "base_value": None, "perturbed_value": None}]
    for param, base in base_values.items():
        for pct in pct_levels:
            for sign in (-1, 1):
                factor = 1.0 + sign * pct / 100.0
                rows.append({
                    "run_id": f"{sample_id}__{param}__{'m' if sign < 0 else 'p'}{pct:g}",
                    "sample_id": sample_id, "param": param, "pct": sign * pct,
                    "factor": round(factor, 4), "base_value": base,
                    "perturbed_value": round(base * factor, 6)})
    return pd.DataFrame(rows)


def representative_samples(tumors: dict, n: int) -> list[str]:
    """A few tumors spanning the high-grade axis (deterministic), so a sensitivity sweep
    covers both regimes without running the whole cohort (the ~49-runs-each trap)."""
    hi = sorted(s for s, t in tumors.items() if t.get("high_grade_regime"))
    lo = sorted(s for s, t in tumors.items() if not t.get("high_grade_regime"))
    picks, i = [], 0
    while len(picks) < n and (hi or lo):
        pool = (lo, hi)[i % 2] or (hi, lo)[i % 2]
        if pool:
            picks.append(pool.pop(0))
        i += 1
    return picks


def run_virtual_cohort(cfg, out_dir) -> None:
    """Population sweep: draw correlated synthetic tumors from config/joint_priors.yaml
    (levers + couplings), derive extrinsic axes via the measured transfer, map to PhysiCell params.
    Distinct from the OAT per-config sensitivity (PLAN §1). Writes results/abm/virtual_cohort/draws.csv."""
    jp_path = resolve_path(cfg, "config/joint_priors.yaml")
    ts_path = resolve_path(cfg, "results/couplings/tumor_scores.csv")
    corr_path = resolve_path(cfg, "results/couplings/network_tumorB_marginal_matrix.csv")
    for p in (jp_path, ts_path, corr_path):
        if not p.exists():
            raise SystemExit(f"[virtual_cohort] missing {p} — run WS1 (coupling_core.R + "
                             "19_bifurcation_transfer.R) first")
    jp = yaml.safe_load(jp_path.read_text())
    scores = pd.read_csv(ts_path)
    corr = pd.read_csv(corr_path, index_col=0)
    levers = list(jp["levers"].keys())
    n_draws = int(jp.get("sweep", {}).get("n_draws", 256))
    seed = int(cfg["phase_c"].get("seed", jp["provenance"]["seed"]))

    draws = generate_virtual_cohort(scores, corr, levers, n_draws, seed)
    o2p = cfg["phase_c"].get("omics_to_params", {})
    full = apply_transfer(draws, jp.get("transfer", {}), o2p, dict(cfg["phase_c"]["uq"]["params"]))
    d = ensure_dir(out_dir / "virtual_cohort")
    full.to_csv(d / "draws.csv", index=False)
    print(f"[ok] virtual cohort: {n_draws} correlated synthetic tumors -> {d/'draws.csv'}")
    print(f"[info] levers={levers}")
    print("[info] param spread (min..max): "
          f"prolif_mult {full.proliferation_mult.min():.2f}..{full.proliferation_mult.max():.2f}, "
          f"pressure_hm {full.pressure_half_max.min():.2f}..{full.pressure_half_max.max():.2f}, "
          f"adhesion {full.adhesion_mult.min():.2f}..{full.adhesion_mult.max():.2f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["oat", "virtual_cohort"], default="oat",
                    help="oat = per-config sensitivity (default); virtual_cohort = population draw")
    ap.add_argument("--sample", nargs="*", default=None, help="explicit sample_id(s) (oat)")
    ap.add_argument("--all", action="store_true", help="sweep every tumor (expensive, oat)")
    args = ap.parse_args()

    setup_logging()
    cfg = load_config()
    out_dir = resolve_path(cfg, "results/abm")
    if args.mode == "virtual_cohort":
        run_virtual_cohort(cfg, out_dir)
        return

    uq = cfg["phase_c"]["uq"]
    pct_levels = list(uq["pct_levels"])
    base_params = dict(uq["params"])

    abm = yaml.safe_load(
        resolve_path(cfg, "results/abm/positives_to_physicell.yaml").read_text())
    tumors = abm.get("tumors", {})
    if args.sample:
        samples = args.sample
    elif args.all:
        samples = list(tumors)
    else:
        samples = representative_samples(tumors, int(uq.get("n_representative", 3)))
        print(f"[info] representative UQ on {samples} (use --all to sweep every tumor)")

    total = 0
    for sid in samples:
        tumor = abm["tumors"].get(sid)
        if tumor is None:
            continue
        vals = dict(base_params)
        # per-tumor adhesion base if available (blastemal as representative)
        adh = tumor["cell_types"].get("blastemal", {}).get("adhesion_strength")
        if adh is not None:
            vals["adhesion_strength"] = float(adh)
        df = generate_sweep(sid, vals, pct_levels)
        d = ensure_dir(out_dir / sid / "uq")
        df.to_csv(d / "sweep_manifest.csv", index=False)
        total += len(df)
    n_per = 1 + len(base_params) * len(pct_levels) * 2
    print(f"[ok] UQ manifests for {len(samples)} tumors "
          f"({n_per} runs/tumor incl. replicates-to-be, {total} rows total)")
    print(f"[info] pct levels {pct_levels} x {len(base_params)} params (OAT + baseline)")


if __name__ == "__main__":
    main()
