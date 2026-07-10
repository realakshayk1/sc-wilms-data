#!/usr/bin/env python3
"""ABM-1: map the rigor positives -> per-tumor PhysiCell parameters.

Turns the validated signals into concrete agent-based-model inputs, per tumor:

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

Phase C extension (see config/abm_programs.yaml + phase_c.yaml omics_to_params):
  EMT program (epithelial vs mesenchymal)      -> adhesion_strength + migration_speed
    (CDH1/EPCAM vs VIM/CDH2/SNAI1...)              (reciprocal, per-tumor bounded)
  Contact-inhibition program (crowding)        -> pressure->cycle-entry rule HALF-MAX
    (CDKN1A/B, Hippo, mechano; CCND1/MKI67 neg)   (higher -> brake at lower pressure)
  Hypoxia-tolerance program (HIF1A/VEGFA...)    -> oxygen->necrosis rule HALF-MAX
                                                  (higher tolerance -> necrosis at lower O2)
  IGF program (IGF2/IGF1R; 11p15 LOI)  [v1.1]  -> per-cell IGF2 substrate UPTAKE_RATE
                                                  (Wilms analog of the CRPC androgen axis)

Half-maxes are targeted because they dominate ABM QoIs (Johnson et al. Cell 2025, Fig 2E).
Transforms are transparent and bounded (no fitting): mult = clip(1 +/- k*z, lo, hi) on the
per-tumor z-scored program score. New program-score columns are OPTIONAL: absent -> neutral
(multiplier 1.0 / base half-max), so this runs today and lights up when the scores land.
Reads per_tumor_scores.csv (Phase A) + Phase B histology probabilities; writes
results/abm/positives_to_physicell.yaml + per_tumor_abm_params.csv. Honest scope: this
parameterizes initial conditions + response shapes from cross-sectional data; NOT a fitted
dynamical model.
"""
from __future__ import annotations

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

# optional per-tumor program-score columns (config/abm_programs.yaml); absent -> neutral
EMT_EPI_COL, EMT_MES_COL = "emt_epithelial_score", "emt_mesenchymal_score"
CONTACT_COL, HYPOXIA_COL = "contact_inhibition_score", "hypoxia_score"
IGF_COL = "igf_score"                               # v1.1 IGF2-uptake axis


def clip_lin(z, k, lo, hi):
    return float(np.clip(1.0 + k * (0.0 if pd.isna(z) else z), lo, hi))


def zscore_col(df: pd.DataFrame, col: str) -> pd.Series:
    """Cohort z-score of an optional column; all-zeros (neutral) if absent/degenerate."""
    if col not in df.columns:
        return pd.Series(0.0, index=df.index)
    x = pd.to_numeric(df[col], errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=df.index)
    return ((x - x.mean()) / sd).fillna(0.0)


def main():
    setup_logging()
    cfg = load_config()
    root = Path(cfg["root"])

    pt = pd.read_csv(resolve_path(cfg, "results/mechanotypes/per_tumor_scores.csv"))
    physicell = yaml.safe_load((root / "config" / "physicell.yaml").read_text())
    base = physicell["cell_types"]
    phase_c = yaml.safe_load((root / "config" / "phase_c.yaml").read_text())
    o2p = phase_c.get("omics_to_params", {})
    base_hm = phase_c["uq"]["params"]
    base_pressure_hm = float(base_hm["pressure_half_max"])
    base_necrosis_hm = float(base_hm["oxygen_necrosis_half_max"])
    base_cycle_hm = float(base_hm["oxygen_cycle_half_max"])

    # omics->params scalers (bounded, documented) with safe defaults
    emt_adh_k = float(o2p.get("emt_adhesion_k", 0.25))
    emt_mot_k = float(o2p.get("emt_motility_k", 0.40))
    adh_lo, adh_hi = o2p.get("adhesion_bounds", [0.4, 1.6])
    mot_lo, mot_hi = o2p.get("motility_bounds", [0.4, 2.0])
    press_k = float(o2p.get("pressure_halfmax_k", 0.30))
    necr_k = float(o2p.get("necrosis_halfmax_k", 0.30))
    hm_lo, hm_hi = o2p.get("halfmax_bounds", [0.5, 1.5])
    igf_k = float(o2p.get("igf_uptake_k", 0.40))
    up_lo, up_hi = o2p.get("uptake_bounds", [0.4, 2.5])

    # precompute cohort z-scores for the optional program columns (neutral if absent)
    z_epi, z_mes = zscore_col(pt, EMT_EPI_COL), zscore_col(pt, EMT_MES_COL)
    z_contact, z_hypoxia = zscore_col(pt, CONTACT_COL), zscore_col(pt, HYPOXIA_COL)
    z_igf = zscore_col(pt, IGF_COL)
    has_emt = EMT_EPI_COL in pt.columns or EMT_MES_COL in pt.columns
    has_crowd = CONTACT_COL in pt.columns
    has_hypox = HYPOXIA_COL in pt.columns
    has_igf = IGF_COL in pt.columns

    # F5: compartment-resolved EMT (cell x tumor type). Each compartment uses its own
    # <program>_score__<compartment> column when present (17_abm_program_scores.R), else the
    # tumor-level column. So blastemal / epithelial / stromal cells get distinct adhesion+
    # motility from compartment-specific DE, not one shared tumor value.
    def _z_ct(base_col, tumor_series):
        return {ct: (zscore_col(pt, f"{base_col}__{ct}") if f"{base_col}__{ct}" in pt.columns
                     else tumor_series) for ct in base}
    z_epi_ct, z_mes_ct = _z_ct(EMT_EPI_COL, z_epi), _z_ct(EMT_MES_COL, z_mes)
    emt_resolved = any(f"{EMT_MES_COL}__{ct}" in pt.columns for ct in base)

    # Phase B anaplasia probability per tumor (held-out MIL preds); fall back to subdiagnosis
    pred_path = resolve_path(cfg, "results/classifier/phase_b_mil_phikon-v2.predictions.csv")
    anap_prob = {}
    if pred_path.exists():
        pr = pd.read_csv(pred_path)
        anap_prob = dict(zip(pr["sample_id"], pr["p_mil"]))

    rows, abm = [], {"_meta": {
        "description": "Per-tumor PhysiCell parameters derived from rigor positives",
        "scalers": {"prolif_k": PROLIF_K, "apop_k": APOP_K,
                    "highgrade_prolif_bump": HIGHGRADE_PROLIF_BUMP,
                    "highgrade_adhesion_mult": HIGHGRADE_ADHESION_MULT,
                    "emt_adhesion_k": emt_adh_k, "emt_motility_k": emt_mot_k,
                    "pressure_halfmax_k": press_k, "necrosis_halfmax_k": necr_k,
                    "igf_uptake_k": igf_k},
        "program_scores_present": {"emt": bool(has_emt), "crowding": bool(has_crowd),
                                   "hypoxia": bool(has_hypox), "igf": bool(has_igf)},
        "emt_compartment_resolved": bool(emt_resolved),
        "base_rates": base, "base_half_max": {
            "pressure_cycle": base_pressure_hm, "oxygen_necrosis": base_necrosis_hm,
            "oxygen_cycle": base_cycle_hm}}, "tumors": {}}

    for i, r in pt.iterrows():
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
        prolif_extra = HIGHGRADE_PROLIF_BUMP if high_grade else 1.0

        # --- Phase C extrinsic axes (neutral when program scores absent) ---
        # EMT index: mesenchymal minus epithelial -> reciprocal adhesion / motility
        emt_index = float(z_mes.loc[i] - z_epi.loc[i])
        adhesion_mult_emt = float(np.clip(1.0 - emt_adh_k * emt_index, adh_lo, adh_hi))
        motility_mult = float(np.clip(1.0 + emt_mot_k * emt_index, mot_lo, mot_hi))
        adhesion_mult = (HIGHGRADE_ADHESION_MULT if high_grade else 1.0) * adhesion_mult_emt
        # crowding / hypoxia -> rule half-maxes (higher program -> lower half-max)
        pressure_half_max = round(base_pressure_hm
                                  * float(np.clip(1.0 - press_k * z_contact.loc[i], hm_lo, hm_hi)), 4)
        necrosis_half_max = round(base_necrosis_hm
                                  * float(np.clip(1.0 - necr_k * z_hypoxia.loc[i], hm_lo, hm_hi)), 4)
        # IGF2 uptake (v1.1): per-tumor IGF program scales each cell's uptake rate
        igf_uptake_mult = float(np.clip(1.0 + igf_k * z_igf.loc[i], up_lo, up_hi))

        cells = {}
        for ct, b in base.items():
            # compartment-specific EMT -> this cell type's adhesion + motility
            emt_ct = float(z_mes_ct[ct].loc[i] - z_epi_ct[ct].loc[i])
            adh_ct = ((HIGHGRADE_ADHESION_MULT if high_grade else 1.0)
                      * float(np.clip(1.0 - emt_adh_k * emt_ct, adh_lo, adh_hi)))
            mot_ct = float(np.clip(1.0 + emt_mot_k * emt_ct, mot_lo, mot_hi))
            cells[ct] = {
                "proliferation_rate": round(b["proliferation_rate"] * prolif_mult * prolif_extra, 6),
                "apoptosis_rate": round(b["apoptosis_rate"] * apop_mult, 6),
                "adhesion_strength": round(b.get("adhesion_strength", 0.5) * adh_ct, 4),
                "migration_speed": round(b.get("migration_speed", 0.3) * mot_ct, 4),
                "igf_uptake_rate": round(b.get("igf_uptake_rate", 0.001) * igf_uptake_mult, 6),
                "ecm_secretion_rate": round(b.get("ecm_secretion_rate", 0.0), 6),
            }
            if "ecm_stiffness" in b:
                cells[ct]["ecm_stiffness"] = b["ecm_stiffness"]
        abm["tumors"][sid] = {
            "initial_fractions": init_frac, "high_grade_regime": high_grade,
            "half_max": {"pressure_cycle": pressure_half_max,
                         "oxygen_necrosis": necrosis_half_max,
                         "oxygen_cycle": round(base_cycle_hm, 4)},
            "cell_types": cells}
        rows.append({"sample_id": sid, **{f"init_{k}": v for k, v in init_frac.items()},
                     "proliferation_mult": prolif_mult, "apoptosis_mult": apop_mult,
                     "emt_index_z": round(emt_index, 3), "adhesion_mult": round(adhesion_mult, 3),
                     "motility_mult": round(motility_mult, 3),
                     "igf_uptake_mult": round(igf_uptake_mult, 3),
                     "pressure_half_max": pressure_half_max, "necrosis_half_max": necrosis_half_max,
                     "high_grade": high_grade,
                     "anaplasia_prob": None if pd.isna(p_anap) else round(float(p_anap), 3),
                     "relapse": r.get("relapse")})

    out_dir = resolve_path(cfg, "results/abm"); ensure_dir(out_dir)
    (out_dir / "positives_to_physicell.yaml").write_text(yaml.safe_dump(abm, sort_keys=False))
    tab = pd.DataFrame(rows)
    tab.to_csv(out_dir / "per_tumor_abm_params.csv", index=False)
    print(f"[ok] {len(rows)} tumors mapped -> {out_dir/'positives_to_physicell.yaml'}")
    print(f"[ok] table -> {out_dir/'per_tumor_abm_params.csv'}")
    print(f"[info] program scores present: emt={has_emt} crowding={has_crowd} "
          f"hypoxia={has_hypox} igf={has_igf}  (absent -> neutral)")
    # quick sanity: do high-grade / relapse tumors get higher proliferation multipliers?
    if "relapse" in tab and tab["relapse"].notna().any():
        rel = tab.dropna(subset=["relapse"])
        print(f"[check] mean proliferation_mult  relapse={rel[rel.relapse==1]['proliferation_mult'].mean():.2f}  "
              f"no-relapse={rel[rel.relapse==0]['proliferation_mult'].mean():.2f}")
    print(tab.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
