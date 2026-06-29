#!/usr/bin/env python3
"""Negative controls: shuffled labels and random morphology features."""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from classifier_utils import (
    FEATURE_COLS,
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
    seed = set_seed_logged(cfg["features"]["seed"], "negative_controls")
    pb = cfg["phase_b"]
    out_json = resolve_path(cfg, cfg["paths"]["phase_b"]["negative_controls_json"])
    ensure_dir(out_json.parent)

    if out_json.exists() and not args.force:
        print(f"[skip] Negative controls exist: {out_json}")
        return

    df = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"]))
    sig = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"]))
    deconv = sig[["spot_id", "dominant_state"] + [f"deconv_{s}" for s in ["blastemal", "epithelial", "stromal"]]]

    samples = sorted(df["sample_id"].astype(str).unique())
    if len(samples) < 2:
        raise SystemExit("Need >= 2 samples for negative controls")

    holdout = samples[-1]
    train = df[df["sample_id"].astype(str) != holdout].copy()
    test = df[df["sample_id"].astype(str) == holdout].copy()
    rng = np.random.default_rng(seed)

    results: dict = {"holdout_sample": holdout, "seed": seed}

    # Control 1: shuffled weak labels within training set
    train_shuf = train.copy()
    train_shuf["weak_label"] = rng.permutation(train_shuf["weak_label"].values)
    y_pred, _ = train_and_predict(train_shuf, test, pb.get("classifier", "random_forest"), seed)
    test_shuf = test.copy()
    test_shuf["pred_state"] = y_pred
    frac_shuf = spot_fractions_from_preds(test_shuf)
    eval_shuf = evaluate_spot_fractions(frac_shuf, deconv)
    results["shuffled_labels"] = {
        "dominant_state_agreement": eval_shuf["dominant_state_agreement"],
        "correlations": eval_shuf["correlations"],
    }

    # Control 2: random morphology features
    train_rand = train.copy()
    test_rand = test.copy()
    for col in FEATURE_COLS:
        mu, sigma = train[col].mean(), train[col].std() or 1.0
        train_rand[col] = rng.normal(mu, sigma, len(train_rand))
        test_rand[col] = rng.normal(mu, sigma, len(test_rand))
    y_pred_rand, _ = train_and_predict(train_rand, test_rand, pb.get("classifier", "random_forest"), seed + 1)
    test_rand = test_rand.copy()
    test_rand["pred_state"] = y_pred_rand
    frac_rand = spot_fractions_from_preds(test_rand)
    eval_rand = evaluate_spot_fractions(frac_rand, deconv)
    results["random_features"] = {
        "dominant_state_agreement": eval_rand["dominant_state_agreement"],
        "correlations": eval_rand["correlations"],
    }

    # Control 3: real labels (reference on same holdout)
    y_pred_real, _ = train_and_predict(train, test, pb.get("classifier", "random_forest"), seed + 2)
    test_real = test.copy()
    test_real["pred_state"] = y_pred_real
    frac_real = spot_fractions_from_preds(test_real)
    eval_real = evaluate_spot_fractions(frac_real, deconv)
    results["real_labels_reference"] = {
        "dominant_state_agreement": eval_real["dominant_state_agreement"],
        "correlations": eval_real["correlations"],
    }

    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

    print("[ok] Negative controls:")
    print(f"     real epi r={eval_real['correlations']['epithelial']['pearson_r']:.3f}")
    print(f"     shuf epi r={eval_shuf['correlations']['epithelial']['pearson_r']:.3f}")
    print(f"     rand epi r={eval_rand['correlations']['epithelial']['pearson_r']:.3f}")
    print(f"[ok] -> {out_json}")


if __name__ == "__main__":
    main()
