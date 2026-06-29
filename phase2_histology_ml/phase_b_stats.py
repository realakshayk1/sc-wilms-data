#!/usr/bin/env python3
"""Shared validation statistics for Phase B (V-1 hardening).

DeLong variance-based CI for a single ROC AUC (fast, parametric), plus a
label-permutation test on a fixed prediction vector (non-parametric, exact-ish).
Both operate on held-out predictions so they layer onto any LOTO classifier
without retraining.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def _compute_midrank(x: np.ndarray) -> np.ndarray:
    """Midranks (ties averaged) — the kernel of the fast DeLong estimator."""
    order = np.argsort(x)
    x_sorted = x[order]
    n = len(x)
    T = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and x_sorted[j] == x_sorted[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(n, dtype=float)
    T2[order] = T
    return T2


def delong_auc_ci(y_true: np.ndarray, y_score: np.ndarray, alpha: float = 0.05):
    """AUC with DeLong asymptotic variance and (1-alpha) CI.

    Returns dict(auc, var, se, ci_low, ci_high, n_pos, n_neg). CI is computed on
    the logit scale then back-transformed so it stays inside [0, 1].
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    m, n = len(pos), len(neg)
    if m == 0 or n == 0:
        return {"auc": float("nan"), "var": float("nan"), "se": float("nan"),
                "ci_low": float("nan"), "ci_high": float("nan"), "n_pos": m, "n_neg": n}
    tx = _compute_midrank(pos)
    ty = _compute_midrank(neg)
    txy = _compute_midrank(np.concatenate([pos, neg]))
    txy_pos = txy[:m]
    txy_neg = txy[m:]
    auc = (np.sum(txy_pos) - m * (m + 1) / 2.0) / (m * n)
    v01 = (txy_pos - tx) / n            # structural components over positives
    v10 = 1.0 - (txy_neg - ty) / m      # over negatives
    s01 = np.var(v01, ddof=1) if m > 1 else 0.0
    s10 = np.var(v10, ddof=1) if n > 1 else 0.0
    var = s01 / m + s10 / n
    se = float(np.sqrt(var)) if var > 0 else 0.0
    z = stats.norm.ppf(1 - alpha / 2)
    if se > 0 and 0 < auc < 1:
        lg = np.log(auc / (1 - auc))
        lse = se / (auc * (1 - auc))
        lo = 1 / (1 + np.exp(-(lg - z * lse)))
        hi = 1 / (1 + np.exp(-(lg + z * lse)))
    else:
        lo, hi = max(0.0, auc - z * se), min(1.0, auc + z * se)
    return {"auc": float(auc), "var": float(var), "se": se,
            "ci_low": float(lo), "ci_high": float(hi), "n_pos": int(m), "n_neg": int(n)}


def permutation_auc_p(y_true: np.ndarray, y_score: np.ndarray, n_perm: int = 10000,
                      seed: int = 42):
    """One-sided permutation p-value for AUC > 0.5 by shuffling labels.

    Operates on the fixed held-out prediction vector (no model retraining), so it
    tests whether the ranking carries label information. p = (1 + #{auc_perm >=
    auc_obs}) / (n_perm + 1).
    """
    from sklearn.metrics import roc_auc_score
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if len(np.unique(y_true)) < 2:
        return {"auc": float("nan"), "p_value": float("nan"), "n_perm": n_perm}
    auc_obs = roc_auc_score(y_true, y_score)
    rng = np.random.default_rng(seed)
    ge = 0
    yt = y_true.copy()
    for _ in range(n_perm):
        rng.shuffle(yt)
        if roc_auc_score(yt, y_score) >= auc_obs:
            ge += 1
    return {"auc": float(auc_obs), "p_value": (1 + ge) / (n_perm + 1), "n_perm": int(n_perm)}


def delong_paired_test(y_true: np.ndarray, score_a: np.ndarray, score_b: np.ndarray):
    """Paired DeLong test that AUC(A) != AUC(B) on the same samples.

    Returns dict(auc_a, auc_b, delta, z, p_value). Uses the full DeLong covariance
    so it accounts for the correlation between the two score vectors.
    """
    y_true = np.asarray(y_true).astype(int)
    order = np.argsort(-y_true)  # positives first
    y_true = y_true[order]
    scores = np.vstack([np.asarray(score_a, float)[order], np.asarray(score_b, float)[order]])
    m = int(np.sum(y_true == 1))
    n = len(y_true) - m
    if m == 0 or n == 0:
        return {"auc_a": float("nan"), "auc_b": float("nan"), "delta": float("nan"),
                "z": float("nan"), "p_value": float("nan")}
    pos = scores[:, :m]
    neg = scores[:, m:]
    k = 2
    tx = np.array([_compute_midrank(pos[r]) for r in range(k)])
    ty = np.array([_compute_midrank(neg[r]) for r in range(k)])
    tz = np.array([_compute_midrank(scores[r]) for r in range(k)])
    aucs = (np.sum(tz[:, :m], axis=1) - m * (m + 1) / 2.0) / (m * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    cov = sx / m + sy / n
    cov = np.atleast_2d(cov)
    delta = aucs[0] - aucs[1]
    L = np.array([1.0, -1.0])
    var = L @ cov @ L
    se = np.sqrt(var) if var > 0 else 0.0
    z = delta / se if se > 0 else 0.0
    p = 2 * stats.norm.sf(abs(z))
    return {"auc_a": float(aucs[0]), "auc_b": float(aucs[1]), "delta": float(delta),
            "z": float(z), "p_value": float(p)}
