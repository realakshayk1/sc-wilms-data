#!/usr/bin/env python3
"""Publication-style figures for Phase B histology ML."""

from __future__ import annotations

import json
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr

from spatial_utils import CELL_STATES
from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

STATE_COLORS = {
    "blastemal": "#E64B35",
    "epithelial": "#4DBBD5",
    "stromal": "#00A087",
}


def plot_deconv_scatter(merged: pd.DataFrame, fig_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for ax, state in zip(axes, CELL_STATES, strict=True):
        x = merged[f"deconv_{state}"]
        y = merged[f"frac_{state}"]
        r, p = pearsonr(x, y)
        ax.scatter(x, y, alpha=0.35, s=18, c=STATE_COLORS[state], edgecolors="none")
        lims = [0, max(x.max(), y.max()) * 1.05]
        ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f"RNA deconv ({state})")
        ax.set_ylabel(f"H&E fraction ({state})")
        ax.set_title(f"{state.capitalize()}\nPearson r = {r:.2f}, p = {p:.1e}")
        ax.set_aspect("equal")
    fig.suptitle("Visium spot fractions: morphology classifier vs transcriptomic programs", y=1.02, fontsize=12)
    out = fig_dir / "phase_b_deconv_validation.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_dominant_agreement(merged: pd.DataFrame, fig_dir: Path) -> None:
    merged = merged.copy()
    merged["pred_dominant"] = merged[[f"frac_{s}" for s in CELL_STATES]].idxmax(axis=1)
    merged["pred_dominant"] = merged["pred_dominant"].str.replace("frac_", "")
    ct = pd.crosstab(merged["dominant_state"], merged["pred_dominant"])
    for s in CELL_STATES:
        if s not in ct.index:
            ct.loc[s] = 0
        if s not in ct.columns:
            ct[s] = 0
    ct = ct.loc[CELL_STATES, CELL_STATES]

    fig, ax = plt.subplots(figsize=(5, 4.2))
    sns.heatmap(ct, annot=True, fmt="d", cmap="Blues", cbar=False, ax=ax, linewidths=0.5)
    ax.set_xlabel("H&E dominant state (nucleus classifier)")
    ax.set_ylabel("RNA dominant state (spot programs)")
    ax.set_title("Dominant compartment agreement per spot")
    out = fig_dir / "phase_b_dominant_state_confusion.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_fraction_by_histology(frac_df: pd.DataFrame, fig_dir: Path) -> None:
    df = frac_df.copy()
    df["histology"] = df["subdiagnosis"].str.capitalize()
    long = df.melt(
        id_vars=["histology"],
        value_vars=[f"frac_{s}" for s in CELL_STATES],
        var_name="state",
        value_name="fraction",
    )
    long["state"] = long["state"].str.replace("frac_", "")

    fig, ax = plt.subplots(figsize=(7, 4))
    sns.boxplot(data=long, x="state", y="fraction", hue="histology", palette={"Favorable": "#4DBBD5", "Anaplastic": "#E64B35"}, ax=ax)
    ax.set_xlabel("Predicted compartment")
    ax.set_ylabel("Mean spot fraction")
    ax.set_title("H&E-derived composition by tumor histology")
    ax.legend(title="Histology", frameon=False)
    out = fig_dir / "phase_b_fractions_by_histology.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_segmentation_mosaic(cfg: dict, fig_dir: Path, n_tiles: int = 12, seed: int = 42) -> None:
    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    overlay_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["nuclei_dir"]) / "overlays"
    manifest_path = tiles_dir / "tiles_manifest.json"
    if not manifest_path.exists() or not overlay_dir.exists():
        return

    import json as json_mod

    with open(manifest_path) as f:
        manifest = json_mod.load(f)
    rng = np.random.default_rng(seed)
    picks = rng.choice(manifest, size=min(n_tiles, len(manifest)), replace=False)

    cols = 4
    rows = int(np.ceil(len(picks) / cols))
    fig = plt.figure(figsize=(cols * 2.2, rows * 2.2))
    gs = GridSpec(rows, cols, figure=fig, wspace=0.05, hspace=0.12)

    for i, entry in enumerate(picks):
        spot_id = entry["spot_id"]
        overlay_path = overlay_dir / f"{spot_id}_overlay.png"
        if not overlay_path.exists():
            continue
        bgr = cv2.imread(str(overlay_path))
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        ax = fig.add_subplot(gs[i // cols, i % cols])
        ax.imshow(rgb)
        ax.set_title(entry.get("dominant_state", "")[:3], fontsize=8)
        ax.axis("off")

    fig.suptitle("Nuclei segmentation QC (sample Visium spot tiles)", fontsize=11, y=0.98)
    out = fig_dir / "phase_b_segmentation_mosaic.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_classifier_summary(metrics_path: Path, deconv_path: Path, fig_dir: Path) -> None:
    if not metrics_path.exists() or not deconv_path.exists():
        return
    with open(metrics_path) as f:
        metrics = json.load(f)
    with open(deconv_path) as f:
        deconv = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

    # Left: balanced accuracy metrics
    labels = ["Nucleus\nbalanced acc.", "Spot majority\nacc.", "Dominant state\nagreement"]
    vals = [
        metrics.get("balanced_accuracy_nucleus", 0),
        metrics.get("balanced_accuracy_spot_majority", 0),
        deconv.get("dominant_state_agreement", 0),
    ]
    colors = ["#3C5488", "#8491B4", "#00A087"]
    axes[0].bar(labels, vals, color=colors, edgecolor="white")
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Score")
    axes[0].set_title("Classifier performance summary")
    for i, v in enumerate(vals):
        axes[0].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)

    # Right: Pearson r per compartment
    corrs = deconv.get("correlations", {})
    states = CELL_STATES
    r_vals = [corrs.get(s, {}).get("pearson_r", 0) for s in states]
    axes[1].bar(states, r_vals, color=[STATE_COLORS[s] for s in states], edgecolor="white")
    axes[1].set_ylim(0, max(0.6, max(r_vals) * 1.2))
    axes[1].set_ylabel("Pearson r")
    axes[1].set_title("H&E vs RNA fraction correlation")
    for i, v in enumerate(r_vals):
        axes[1].text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)

    fig.tight_layout()
    out = fig_dir / "phase_b_classifier_summary.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_loto_summary(loto_path: Path, fig_dir: Path) -> None:
    if not loto_path.exists():
        return
    with open(loto_path) as f:
        loto = json.load(f)
    folds = loto.get("folds", [])
    if not folds:
        return
    states = CELL_STATES
    mean_r = loto.get("mean_pearson_r", {})
    std_r = loto.get("std_pearson_r", {})
    fig, ax = plt.subplots(figsize=(6, 4))
    vals = [mean_r.get(s, 0) for s in states]
    errs = [std_r.get(s, 0) for s in states]
    ax.bar(states, vals, yerr=errs, color=[STATE_COLORS[s] for s in states], capsize=4, edgecolor="white")
    ax.set_ylabel("Pearson r (mean ± SD across LOTO folds)")
    ax.set_title(f"LOTO validation ({loto.get('n_folds', 0)} folds)")
    ax.set_ylim(0, max(0.6, max(vals) * 1.3))
    out = fig_dir / "phase_b_loto_validation.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_negative_controls(neg_path: Path, fig_dir: Path) -> None:
    if not neg_path.exists():
        return
    with open(neg_path) as f:
        neg = json.load(f)
    labels = ["real_labels_reference", "shuffled_labels", "random_features"]
    epi_r = [neg.get(k, {}).get("correlations", {}).get("epithelial", {}).get("pearson_r", 0) for k in labels]
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.bar(["Real", "Shuffled labels", "Random features"], epi_r, color=["#00A087", "#E64B35", "#B09C85"])
    ax.set_ylabel("Epithelial fraction Pearson r")
    ax.set_title("Negative controls vs real morphology model")
    ax.set_ylim(0, max(0.5, max(epi_r) * 1.2))
    out = fig_dir / "phase_b_negative_controls.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def main() -> None:
    setup_logging()
    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "figures")
    fig_dir = ensure_dir(resolve_path(cfg, cfg["paths"]["dirs"]["figures"]))

    sns.set_theme(style="whitegrid", font_scale=0.95)

    frac_df = pd.read_csv(resolve_path(cfg, cfg["paths"]["phase_b"]["spot_fractions_csv"]))
    sig_df = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"]))
    merged = frac_df.merge(sig_df[["spot_id", "dominant_state"] + [f"deconv_{s}" for s in CELL_STATES]], on="spot_id")

    plot_deconv_scatter(merged, fig_dir)
    plot_dominant_agreement(merged, fig_dir)
    plot_fraction_by_histology(frac_df, fig_dir)
    plot_segmentation_mosaic(cfg, fig_dir, seed=seed)
    plot_classifier_summary(
        resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / "classifier_metrics.json",
        resolve_path(cfg, cfg["paths"]["phase_b"]["deconv_comparison_csv"]),
        fig_dir,
    )
    plot_loto_summary(resolve_path(cfg, cfg["paths"]["phase_b"]["loto_json"]), fig_dir)
    plot_negative_controls(resolve_path(cfg, cfg["paths"]["phase_b"]["negative_controls_json"]), fig_dir)
    print(f"[ok] Phase B figures -> {fig_dir}")


if __name__ == "__main__":
    main()
