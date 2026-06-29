#!/usr/bin/env python3
"""ABM-1: map the rigor positives -> per-tumor PhysiCell parameters.

Turns the four validated signals into concrete agent-based-model inputs, per tumor:

  POSITIVE (evidence)                          -> PhysiCell parameter
  ----------------------------------------------------------------------------------
  Compartment composition (12_composition;     -> initial cell-type FRACTIONS
    epithelial^ anaplastic FDR<0.05)              (blastemal/epithelial/stromal)
  Proliferation/E2F/G2M^ in relapse            -> proliferation_RATE multiplier
    (15_hallmark_gsea padj~1e-29; A-3 OR~4)       (per-tumor, bounded)
  p53-target activity (13/14 moderated DE;      -> apoptosis_RATE multiplier
    TP53 targets v relapse axis)                  (low p53 activity -> less apoptosis)
  Anaplasia from H&E (14/15 Phikon AUC 0.73,    -> high-grade REGIME flag
    16 StarDist morphology 0.69)                  (extra proliferation, less adhesion)

Transforms are transparent and bounded (no fitting): rate_mult = clip(1 + k*z, lo, hi)
on the per-tumor z-scored score. Reads per_tumor_scores.csv (Phase A) + the Phase B
histology probabilities; writes results/abm/positives_to_physicell.yaml +
per_tumor_abm_params.csv. Honest scope: this parameterizes initial conditions from
cross-sectional data; it is NOT a fitted dynamical model.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from utils import ensure_dir, load_config, resolve_path, setup_logging

# bounded linear scalers on per-tumor z-scores (slope k, clip lo/hi) — documented, not fit
PROLIF_K, PROLIF_LO, PROLIF_HI = 0.60, 0.40, 2.50   # +1 SD proliferation -> 1.6x rate
APOP_K, APOP_LO, APOP_HI = 0.50, 0.40, 2.00         # +1 SD p53 activity  -> 1.5x apoptosis
HIGHGRADE_PROLIF_BUMP = 1.25                        # anaplastic regime extra proliferation
HIGHGRADE_ADHESION_MULT = 0.80                      # anaplastic regime reduced adhesion


def clip_lin(z, k, lo, hi):
    return float(np.clip(1.0 + k * (0.0 if pd.isna(z) else z), lo, hi))


def main():
    setup_logging()
    cfg = load_config()
    root = Path(cfg["root"])

    pt = pd.read_csv(resolve_path(cfg, "results/mechanotypes/per_tumor_scores.csv"))
    physicell = yaml.safe_load((root / "config" / "physicell.yaml").read_text())
    base = physicell["cell_types"]

    # Phase B anaplasia probability per tumor (held-out MIL preds); fall back to subdiagnosis
    pred_path = resolve_path(cfg, "results/classifier/phase_b_mil_phikon-v2.predictions.csv")
    anap_prob = {}
    if pred_path.exists():
        pr = pd.read_csv(pred_path)
        anap_prob = dict(zip(pr["sample_id"], pr["p_mil"]))

    comp_cols = [c for c in pt.columns if c.endswith("_frac")]
    rows, abm = [], {"_meta": {
        "description": "Per-tumor PhysiCell parameters derived from rigor positives",
        "scalers": {"prolif_k": PROLIF_K, "apop_k": APOP_K,
                    "highgrade_prolif_bump": HIGHGRADE_PROLIF_BUMP,
                    "highgrade_adhesion_mult": HIGHGRADE_ADHESION_MULT},
        "base_rates": base}, "tumors": {}}

    for _, r in pt.iterrows():
        sid = r["sample_id"]
        # initial fractions (renormalized over the three compartments)
        fr = np.array([r.get(f"{c}_frac", np.nan) for c in ["blastemal", "epithelial", "stromal"]], float)
        if np.all(np.isnan(fr)):
            continue
        fr = np.nan_to_num(fr, nan=0.0)
        fr = fr / fr.sum() if fr.sum() > 0 else np.array([1 / 3, 1 / 3, 1 / 3])
        init_frac = dict(zip(["blastemal", "epithelial", "stromal"], np.round(fr, 4).tolist()))

        prolif_mult = clip_lin(r.get("proliferation_score"), PROLIF_K, PROLIF_LO, PROLIF_HI)
        # low p53-target activity -> LESS apoptosis: scale around the score (high score -> more)
        apop_mult = clip_lin(r.get("tp53_target_score"), APOP_K, APOP_LO, APOP_HI)

        p_anap = anap_prob.get(sid, np.nan)
        sub = str(r.get("subdiagnosis", ""))
        high_grade = bool((p_anap >= 0.5) if not pd.isna(p_anap) else (sub == "anaplastic"))
        adhesion_mult = HIGHGRADE_ADHESION_MULT if high_grade else 1.0
        prolif_extra = HIGHGRADE_PROLIF_BUMP if high_grade else 1.0

        cells = {}
        for ct, b in base.items():
            cells[ct] = {
                "proliferation_rate": round(b["proliferation_rate"] * prolif_mult * prolif_extra, 6),
                "apoptosis_rate": round(b["apoptosis_rate"] * apop_mult, 6),
                "adhesion_strength": round(b.get("adhesion_strength", 0.5) * adhesion_mult, 4),
            }
            if "ecm_stiffness" in b:
                cells[ct]["ecm_stiffness"] = b["ecm_stiffness"]
        abm["tumors"][sid] = {"initial_fractions": init_frac, "high_grade_regime": high_grade,
                              "cell_types": cells}
        rows.append({"sample_id": sid, **{f"init_{k}": v for k, v in init_frac.items()},
                     "proliferation_mult": prolif_mult, "apoptosis_mult": apop_mult,
                     "high_grade": high_grade, "anaplasia_prob": None if pd.isna(p_anap) else round(float(p_anap), 3),
                     "relapse": r.get("relapse")})

    out_dir = resolve_path(cfg, "results/abm"); ensure_dir(out_dir)
    (out_dir / "positives_to_physicell.yaml").write_text(yaml.safe_dump(abm, sort_keys=False))
    tab = pd.DataFrame(rows)
    tab.to_csv(out_dir / "per_tumor_abm_params.csv", index=False)
    print(f"[ok] {len(rows)} tumors mapped -> {out_dir/'positives_to_physicell.yaml'}")
    print(f"[ok] table -> {out_dir/'per_tumor_abm_params.csv'}")
    # quick sanity: do high-grade / relapse tumors get higher proliferation multipliers?
    if "relapse" in tab and tab["relapse"].notna().any():
        rel = tab.dropna(subset=["relapse"])
        print(f"[check] mean proliferation_mult  relapse={rel[rel.relapse==1]['proliferation_mult'].mean():.2f}  "
              f"no-relapse={rel[rel.relapse==0]['proliferation_mult'].mean():.2f}")
    print(tab.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
