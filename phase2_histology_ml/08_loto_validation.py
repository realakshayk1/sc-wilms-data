#!/usr/bin/env python3
"""Leave-one-tumor-out (LOTO) validation for Phase B morphology classifier."""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from classifier_utils import (
    confident_spot_ids,
    evaluate_spot_fractions,
    spot_fractions_from_preds,
    train_and_predict,
)
from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "loto_validation")
    pb = cfg["phase_b"]
    val = pb.get("validation", {})
    min_nuclei = int(val.get("loto_min_nuclei", 200))

    features_path = resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"])
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    out_json = resolve_path(cfg, cfg["paths"]["phase_b"]["loto_json"])
    ensure_dir(out_json.parent)

    if out_json.exists() and not args.force:
        print(f"[skip] LOTO results exist: {out_json}")
        return

    df = pd.read_parquet(features_path)
    sig = pd.read_parquet(sig_path)
    confident = confident_spot_ids(sig)
    df = df[df["spot_id"].isin(confident)].copy()

    samples = sorted(df["sample_id"].astype(str).unique())
    if len(samples) < 2:
        raise SystemExit("Need >= 2 samples for LOTO validation")

    deconv = sig[["spot_id", "dominant_state"] + [f"deconv_{s}" for s in ["blastemal", "epithelial", "stromal"]]]
    fold_results = []

    for holdout in samples:
        train = df[df["sample_id"].astype(str) != holdout]
        test = df[df["sample_id"].astype(str) == holdout]
        if len(train) < min_nuclei or len(test) < 30:
            print(f"[skip] {holdout}: train={len(train)} test={len(test)}")
            continue

        y_pred, metrics = train_and_predict(
            train, test, pb.get("classifier", "random_forest"), seed
        )
        test = test.copy()
        test["pred_state"] = y_pred
        frac_df = spot_fractions_from_preds(test)
        eval_out = evaluate_spot_fractions(frac_df, deconv)

        fold_results.append(
            {
                "holdout_sample": holdout,
                "subdiagnosis": str(test["subdiagnosis"].iloc[0]) if "subdiagnosis" in test.columns else "",
                "n_test_nuclei": int(len(test)),
                "n_test_spots": int(test["spot_id"].nunique()),
                "balanced_accuracy_nucleus": metrics["balanced_accuracy_nucleus"],
                "per_class_balanced_accuracy": metrics["per_class"],
                "dominant_state_agreement": eval_out["dominant_state_agreement"],
                "correlations": eval_out["correlations"],
            }
        )
        print(
            f"[loto] {holdout}: dom_agree={eval_out['dominant_state_agreement']:.3f} "
            f"epi_r={eval_out['correlations']['epithelial']['pearson_r']:.3f}"
        )

    if not fold_results:
        raise SystemExit("No LOTO folds completed — check min_nuclei threshold")

    def _mean_metric(key: str, subkey: str | None = None) -> float:
        vals = []
        for f in fold_results:
            if subkey:
                vals.append(f["correlations"][subkey][key])
            else:
                vals.append(f[key])
        return float(np.mean(vals))

    summary = {
        "n_folds": len(fold_results),
        "mean_dominant_state_agreement": _mean_metric("dominant_state_agreement"),
        "mean_balanced_accuracy_nucleus": _mean_metric("balanced_accuracy_nucleus"),
        "mean_pearson_r": {
            state: _mean_metric("pearson_r", state) for state in ["blastemal", "epithelial", "stromal"]
        },
        "std_pearson_r": {
            state: float(np.std([f["correlations"][state]["pearson_r"] for f in fold_results]))
            for state in ["blastemal", "epithelial", "stromal"]
        },
        "folds": fold_results,
        "seed": seed,
        "split": "leave-one-sample-out",
    }
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[ok] LOTO validation -> {out_json} ({len(fold_results)} folds)")


if __name__ == "__main__":
    main()
