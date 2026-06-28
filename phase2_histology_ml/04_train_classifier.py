#!/usr/bin/env python3
"""FR-B5: Train RF on morphology; spot-level holdout by sample; balanced accuracy."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, classification_report
from sklearn.model_selection import GroupShuffleSplit

from spatial_utils import CELL_STATES
from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

FEATURE_COLS = [
    "area",
    "eccentricity",
    "solidity",
    "major_axis_length",
    "texture_var",
    "hematoxylin_intensity",
    "neighbor_density",
]


def confident_spot_ids(sig_path: Path, min_margin: float = 0.12) -> set[str]:
    """Keep spots where dominant program clearly exceeds the runner-up."""
    sig = pd.read_parquet(sig_path)
    cols = [f"deconv_{s}" for s in CELL_STATES]
    vals = sig[cols].to_numpy()
    top2 = np.sort(vals, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]
    keep = sig.loc[margin >= min_margin, "spot_id"]
    return set(keep.astype(str))


def build_classifier(name: str, seed: int):
    if name == "gradient_boosting":
        return GradientBoostingClassifier(random_state=seed)
    return RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
        min_samples_leaf=2,
    )


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "train_classifier")
    pb = cfg["phase_b"]

    in_parquet = resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"])
    model_path = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"])
    metrics_path = model_path.parent / "classifier_metrics.json"
    methods_path = resolve_path(cfg, cfg["paths"]["phase_b"]["methods_json"])
    ensure_dir(model_path.parent)

    if model_path.exists() and not args.force:
        print(f"[skip] Model exists: {model_path}")
        return

    if not in_parquet.exists():
        raise SystemExit("Run 03_nucleus_features.py first")

    df = pd.read_parquet(in_parquet)
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    if sig_path.exists():
        confident = confident_spot_ids(sig_path)
        before = len(df)
        df = df[df["spot_id"].isin(confident)].copy()
        print(f"[filter] high-confidence spots: {len(confident)} -> {len(df)} nuclei (from {before})")

    if len(df) < 30:
        raise SystemExit(f"Too few nuclei for training: {len(df)}")

    # Hold out entire samples (slides) — standard practice to avoid spot leakage
    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=float(pb.get("holdout_sample_fraction", 0.25)),
        random_state=seed,
    )
    groups = df["sample_id"].astype(str)
    train_idx, test_idx = next(splitter.split(df, groups=groups))

    train = df.iloc[train_idx]
    test = df.iloc[test_idx]
    X_train, y_train = train[FEATURE_COLS].values, train["weak_label"].values
    X_test, y_test = test[FEATURE_COLS].values, test["weak_label"].values

    clf = build_classifier(pb.get("classifier", "random_forest"), seed)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    bal_acc = balanced_accuracy_score(y_test, y_pred)
    classes = sorted(np.unique(np.concatenate([y_train, y_test])))
    per_class_dict = {
        c: float(balanced_accuracy_score(y_test == c, y_pred == c))
        for c in classes
    }

    # Spot-level majority-vote accuracy (closer to biological unit of validation)
    test = test.copy()
    test["pred_state"] = y_pred
    spot_acc_rows = []
    for spot_id, grp in test.groupby("spot_id"):
        maj_pred = grp["pred_state"].mode().iloc[0]
        maj_true = grp["weak_label"].mode().iloc[0]
        spot_acc_rows.append(maj_pred == maj_true)
    spot_bal_acc = float(np.mean(spot_acc_rows)) if spot_acc_rows else float("nan")

    metrics = {
        "balanced_accuracy_nucleus": float(bal_acc),
        "balanced_accuracy_spot_majority": spot_bal_acc,
        "per_class_balanced_accuracy": per_class_dict,
        "classification_report": classification_report(y_test, y_pred, output_dict=True),
        "n_train_nuclei": int(len(y_train)),
        "n_test_nuclei": int(len(y_test)),
        "n_train_samples": int(train["sample_id"].nunique()),
        "n_test_samples": int(test["sample_id"].nunique()),
        "holdout_samples": sorted(test["sample_id"].unique().tolist()),
        "seed": seed,
        "split": "GroupShuffleSplit by sample_id",
        "classifier": pb.get("classifier", "random_forest"),
    }

    with open(model_path, "wb") as f:
        pickle.dump({"model": clf, "feature_cols": FEATURE_COLS}, f)
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    methods = {
        "weak_supervision": "Dominant Phase A program per Visium spot; train on high-confidence spots (program margin >= 0.12)",
        "train_test_split": "Hold out entire tumor samples (no spot leakage across slides)",
        "classifier": metrics["classifier"],
        "metrics_primary": "per_class_balanced_accuracy on held-out nuclei",
        "metrics_secondary": "spot-level majority-vote agreement with weak labels",
        "segmentation": pb.get("segmentation_method", "stardist"),
        "stain_normalization": pb.get("stain_normalization", "macenko"),
    }
    with open(methods_path, "w") as f:
        json.dump(methods, f, indent=2)

    print(f"[ok] Nucleus balanced accuracy = {bal_acc:.3f}")
    print(f"[ok] Spot majority accuracy   = {spot_bal_acc:.3f}")
    for c, a in per_class_dict.items():
        print(f"     {c}: {a:.3f}")
    print(f"[ok] Model -> {model_path}")


if __name__ == "__main__":
    main()
