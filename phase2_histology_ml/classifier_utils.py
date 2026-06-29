"""Shared classifier training and spot-level evaluation for Phase B validation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score

from spatial_utils import CELL_STATES

FEATURE_COLS = [
    "area",
    "eccentricity",
    "solidity",
    "major_axis_length",
    "texture_var",
    "hematoxylin_intensity",
    "neighbor_density",
]


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


def confident_spot_ids(sig: pd.DataFrame, min_margin: float = 0.12) -> set[str]:
    cols = [f"deconv_{s}" for s in CELL_STATES]
    vals = sig[cols].to_numpy()
    top2 = np.sort(vals, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]
    keep = sig.loc[margin >= min_margin, "spot_id"]
    return set(keep.astype(str))


def spot_fractions_from_preds(df: pd.DataFrame, pred_col: str = "pred_state") -> pd.DataFrame:
    rows = []
    for spot_id, grp in df.groupby("spot_id"):
        counts = grp[pred_col].value_counts(normalize=True)
        meta = grp.iloc[0]
        row = {
            "spot_id": spot_id,
            "sample_id": meta.get("sample_id", ""),
            "library_id": meta.get("library_id", ""),
            "subdiagnosis": meta.get("subdiagnosis", ""),
        }
        for state in CELL_STATES:
            row[f"frac_{state}"] = float(counts.get(state, 0.0))
        rows.append(row)
    return pd.DataFrame(rows)


def evaluate_spot_fractions(
    frac_df: pd.DataFrame,
    deconv: pd.DataFrame,
    histology_col: str = "subdiagnosis",
) -> dict[str, Any]:
    merged = frac_df.merge(
        deconv[["spot_id", "dominant_state"] + [f"deconv_{s}" for s in CELL_STATES]],
        on="spot_id",
        how="inner",
    )
    correlations: dict[str, dict[str, float]] = {}
    for state in CELL_STATES:
        r, p = pearsonr(merged[f"frac_{state}"], merged[f"deconv_{state}"])
        rho, _ = spearmanr(merged[f"frac_{state}"], merged[f"deconv_{state}"])
        correlations[state] = {
            "pearson_r": float(r),
            "pearson_p": float(p),
            "spearman_rho": float(rho),
        }

    merged["pred_dominant"] = merged[[f"frac_{s}" for s in CELL_STATES]].idxmax(axis=1)
    merged["pred_dominant"] = merged["pred_dominant"].str.replace("frac_", "")
    dom_agree = float((merged["pred_dominant"] == merged["dominant_state"]).mean())

    by_histology: dict[str, Any] = {}
    if histology_col in merged.columns:
        for hist, grp in merged.groupby(merged[histology_col].astype(str)):
            if len(grp) < 10:
                continue
            hcorr = {}
            for state in CELL_STATES:
                r, p = pearsonr(grp[f"frac_{state}"], grp[f"deconv_{state}"])
                hcorr[state] = {"pearson_r": float(r), "pearson_p": float(p), "n_spots": len(grp)}
            by_histology[hist] = hcorr

    return {
        "correlations": correlations,
        "dominant_state_agreement": dom_agree,
        "n_spots": int(len(merged)),
        "by_histology": by_histology,
        "merged": merged,
    }


def train_and_predict(
    train: pd.DataFrame,
    test: pd.DataFrame,
    classifier_name: str,
    seed: int,
) -> tuple[np.ndarray, dict[str, float]]:
    clf = build_classifier(classifier_name, seed)
    X_train = train[FEATURE_COLS].values
    y_train = train["weak_label"].values
    X_test = test[FEATURE_COLS].values
    y_test = test["weak_label"].values
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    bal_acc = float(balanced_accuracy_score(y_test, y_pred))
    classes = sorted(np.unique(np.concatenate([y_train, y_test])))
    per_class = {
        c: float(balanced_accuracy_score(y_test == c, y_pred == c))
        for c in classes
    }
    return y_pred, {"balanced_accuracy_nucleus": bal_acc, "per_class": per_class}
