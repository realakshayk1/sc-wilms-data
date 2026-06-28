#!/usr/bin/env python3
"""FR-B6: Nucleus predictions -> spot fractions; validate vs transcriptomic deconvolution."""

from __future__ import annotations

import argparse
import json
import pickle

import pandas as pd
from scipy.stats import pearsonr, spearmanr

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

CELL_STATES = ["blastemal", "epithelial", "stromal"]


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "spot_fractions")

    features_path = resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"])
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    model_path = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"])
    out_csv = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_fractions_csv"])
    corr_csv = resolve_path(cfg, cfg["paths"]["phase_b"]["deconv_comparison_csv"])
    ensure_dir(out_csv.parent)

    if out_csv.exists() and not args.force:
        print(f"[skip] Fractions exist: {out_csv}")
        return

    df = pd.read_parquet(features_path)
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    clf = bundle["model"]
    cols = bundle["feature_cols"]

    df["pred_state"] = clf.predict(df[cols].values)

    frac_rows = []
    for spot_id, grp in df.groupby("spot_id"):
        counts = grp["pred_state"].value_counts(normalize=True)
        meta = grp.iloc[0]
        row = {
            "spot_id": spot_id,
            "sample_id": meta.get("sample_id", ""),
            "library_id": meta.get("library_id", ""),
            "subdiagnosis": meta.get("subdiagnosis", ""),
        }
        for state in CELL_STATES:
            row[f"frac_{state}"] = float(counts.get(state, 0.0))
        frac_rows.append(row)
    frac_df = pd.DataFrame(frac_rows)
    frac_df.to_csv(out_csv, index=False)

    if not sig_path.exists():
        raise SystemExit(f"Missing spot signatures: {sig_path}")

    deconv = pd.read_parquet(sig_path)[
        ["spot_id"] + [f"deconv_{s}" for s in CELL_STATES] + ["dominant_state"]
    ]
    merged = frac_df.merge(deconv, on="spot_id", how="inner")

    correlations = {}
    for state in CELL_STATES:
        r, p = pearsonr(merged[f"frac_{state}"], merged[f"deconv_{state}"])
        rho, _ = spearmanr(merged[f"frac_{state}"], merged[f"deconv_{state}"])
        correlations[state] = {
            "pearson_r": float(r),
            "pearson_p": float(p),
            "spearman_rho": float(rho),
        }

    # Dominant-state agreement between image classifier and transcriptome
    merged["pred_dominant"] = merged[[f"frac_{s}" for s in CELL_STATES]].idxmax(axis=1)
    merged["pred_dominant"] = merged["pred_dominant"].str.replace("frac_", "")
    dom_agree = float((merged["pred_dominant"] == merged["dominant_state"]).mean())

    out = {
        "correlations": correlations,
        "dominant_state_agreement": dom_agree,
        "n_spots": int(len(merged)),
        "validation_method": "Softmax Phase A program scores per Visium spot (same genes as snRNA-seq)",
        "seed": seed,
    }
    with open(corr_csv, "w") as f:
        json.dump(out, f, indent=2)

    print(f"[ok] Spot fractions -> {out_csv} ({len(frac_df)} spots)")
    print(f"[ok] Dominant-state agreement (H&E vs RNA): {dom_agree:.3f}")
    for state, vals in correlations.items():
        print(f"     {state}: Pearson r={vals['pearson_r']:.3f}, Spearman rho={vals['spearman_rho']:.3f}")


if __name__ == "__main__":
    main()
