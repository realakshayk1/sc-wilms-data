#!/usr/bin/env python3
"""Result figures for Phase A (GSEA + moderated DE), Phase B (histology AUC), and the ABM mapping.

Reads the result artifacts and renders three PNGs into results/figures/:
  phase_a_gsea_de.png  — Hallmark GSEA (relapse axis) + moderated-DE FDR gene counts
  phase_b_histology_auc.png   — tumor-level histology AUC forest with DeLong 95% CIs
  abm_parameters.png    — ABM proliferation multiplier by relapse (biology -> parameter)
All values are read from disk; nothing is hard-coded except reference baselines.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import load_config, resolve_path

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 130, "savefig.bbox": "tight"})
C_UP, C_DN, C_EMB, C_MORPH, C_REF = "#c0392b", "#2471a3", "#1f7a4d", "#b9770e", "#7f8c8d"


def fig_phase_a(cfg, figdir):
    gs = pd.read_csv(resolve_path(cfg, "results/mechanotypes/hallmark_gsea.csv"))
    de = pd.read_csv(resolve_path(cfg, "results/mechanotypes/moderated_de_summary.csv"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.35, 1]})

    # --- (a) relapse-axis Hallmark GSEA, top by FDR ---
    g = gs[(gs.scope == "overall") & (gs.contrast == "relapse_vs_norelapse")].copy()
    g = g.reindex(g["padj"].nsmallest(12).index).iloc[::-1]
    y = np.arange(len(g))
    colors = [C_UP if n > 0 else C_DN for n in g["NES"]]
    ax1.barh(y, g["NES"], color=colors, alpha=0.9)
    ax1.set_yticks(y)
    ax1.set_yticklabels([p.replace("HALLMARK_", "").replace("_", " ").title() for p in g["pathway"]], fontsize=8.5)
    ax1.axvline(0, color="k", lw=0.8)
    for yi, (nes, padj) in enumerate(zip(g["NES"], g["padj"])):
        ax1.text(nes + (0.06 if nes > 0 else -0.06), yi, f"q={padj:.0e}",
                 va="center", ha="left" if nes > 0 else "right", fontsize=7.2, color="#333")
    ax1.set_xlabel("Normalized enrichment score (NES)")
    ax1.set_title("Hallmark GSEA · relapse vs no-relapse\n(↑ = up in relapse; q = BH-FDR)",
                  fontsize=10.5, loc="left")
    ax1.set_xlim(g["NES"].min() - 1.4, g["NES"].max() + 1.4)

    # --- (b) moderated-DE FDR<0.05 gene counts per contrast ---
    de = de.copy()
    de["label"] = de["scope"] + "\n" + de["contrast"].str.replace("_vs_", " v ").str.replace("_", " ")
    de = de.sort_values("edgeR_fdr05", ascending=True)
    y = np.arange(len(de)); h = 0.38
    ax2.barh(y + h / 2, de["edgeR_fdr05"], h, color=C_UP, alpha=0.9, label="edgeR-QLF")
    ax2.barh(y - h / 2, de["voom_fdr05"], h, color=C_DN, alpha=0.7, label="limma-voom")
    ax2.set_yticks(y); ax2.set_yticklabels(de["label"], fontsize=7.8)
    ax2.set_xlabel("genes at FDR < 0.05")
    ax2.set_title("Moderated pseudobulk DE\n(genes at FDR < 0.05)", fontsize=10.5, loc="left")
    ax2.legend(fontsize=8, loc="lower right", frameon=False)
    fig.tight_layout()
    out = figdir / "phase_a_gsea_de.png"; fig.savefig(out); plt.close(fig)
    return out


def fig_phase_b(cfg, figdir):
    mil = json.loads(resolve_path(cfg, "results/classifier/phase_b_mil_phikon-v2.json").read_text())
    sd = json.loads(resolve_path(cfg, "results/classifier/stardist_morphology.json").read_text())
    m = mil["models"]; s = sd["models"]
    # (label, auc, lo, hi, color, is_ref)
    rows = [
        ("Watershed morphology", sd["baseline_watershed_morph_auc_ref"], None, None, C_REF, True),
        ("StarDist morphology", s["stardist_morphology_rf"]["auc"], s["stardist_morphology_rf"]["ci_low"], s["stardist_morphology_rf"]["ci_high"], C_MORPH, False),
        ("Phikon-v1 mean-pool, 60 spots", mil["baseline_v1_meanpool_auc_ref"], None, None, C_REF, True),
        ("Phikon-v2 mean-pool, 200 spots", m["meanpool_logistic"]["auc"], m["meanpool_logistic"]["ci_low"], m["meanpool_logistic"]["ci_high"], C_EMB, False),
        ("Phikon-v2 attention-MIL", m["attention_mil"]["auc"], m["attention_mil"]["ci_low"], m["attention_mil"]["ci_high"], C_EMB, False),
        ("Ensemble: morphology + embedding", s["ensemble_morph_plus_embedding"]["auc"], s["ensemble_morph_plus_embedding"]["ci_low"], s["ensemble_morph_plus_embedding"]["ci_high"], C_UP, False),
    ]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, auc, lo, hi, col, ref) in zip(y, rows):
        if lo is not None:
            ax.plot([lo, hi], [yi, yi], color=col, lw=2.2, alpha=0.85, zorder=2)
        ax.scatter([auc], [yi], s=85, color=col, zorder=3, marker="D" if ref else "o",
                   edgecolor="white", linewidth=0.8)
        ax.text(auc, yi + 0.22, f"{auc:.3f}", ha="center", va="bottom", fontsize=8.5, color=col)
    ax.axvline(0.5, color="k", ls="--", lw=0.9, alpha=0.6)
    ax.text(0.492, 2.5, "chance", fontsize=8, ha="right", va="center", color="k", rotation=90)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("tumor-level histology AUC (leave-one-tumor-out)  ·  bars = DeLong 95% CI")
    ax.set_xlim(0.30, 1.0)
    ax.set_title("Phase B · anaplasia from H&E — real but resolution-capped at ~0.73\n"
                 "MIL vs mean-pool paired DeLong p=0.83 · ensemble vs embedding p=0.57", fontsize=10.5, loc="left")
    fig.tight_layout()
    out = figdir / "phase_b_histology_auc.png"; fig.savefig(out); plt.close(fig)
    return out


def fig_abm(cfg, figdir):
    t = pd.read_csv(resolve_path(cfg, "results/abm/per_tumor_abm_params.csv"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.3]})
    rel = t.dropna(subset=["relapse"])
    groups = [("no-relapse", rel[rel.relapse == 0]["proliferation_mult"], C_DN),
              ("relapse", rel[rel.relapse == 1]["proliferation_mult"], C_UP)]
    for i, (lab, vals, col) in enumerate(groups):
        x = np.random.default_rng(0).normal(i, 0.05, len(vals))
        ax1.scatter(x, vals, color=col, alpha=0.75, s=36, edgecolor="white", linewidth=0.5)
        ax1.plot([i - 0.22, i + 0.22], [vals.mean()] * 2, color=col, lw=2.5)
        ax1.text(i, vals.mean() + 0.04, f"mean {vals.mean():.2f}", ha="center", fontsize=8.5, color=col)
    ax1.axhline(1.0, color=C_REF, ls="--", lw=0.9)
    ax1.set_xticks([0, 1]); ax1.set_xticklabels(["no-relapse", "relapse"])
    ax1.set_ylabel("ABM proliferation-rate multiplier")
    ax1.set_title("ABM proliferation parameter\nby relapse status", fontsize=10.5, loc="left")

    # initial compartment fractions, tumors sorted by epithelial (composition positive)
    ti = t.sort_values("init_epithelial").reset_index(drop=True)
    x = np.arange(len(ti)); bottom = np.zeros(len(ti))
    for comp, col in [("init_blastemal", "#8e44ad"), ("init_epithelial", "#16a085"), ("init_stromal", "#d35400")]:
        ax2.bar(x, ti[comp], bottom=bottom, color=col, width=1.0, label=comp.replace("init_", ""))
        bottom += ti[comp].values
    ax2.set_xlim(-0.5, len(ti) - 0.5); ax2.set_ylim(0, 1)
    ax2.set_xlabel("tumors (sorted by epithelial fraction)"); ax2.set_ylabel("initial cell-type fraction")
    ax2.set_title("Per-tumor initial conditions\nfrom compartment composition", fontsize=10.5, loc="left")
    ax2.legend(fontsize=8, ncol=3, loc="upper center", frameon=False, bbox_to_anchor=(0.5, -0.18))
    fig.tight_layout()
    out = figdir / "abm_parameters.png"; fig.savefig(out); plt.close(fig)
    return out


def main():
    cfg = load_config()
    figdir = resolve_path(cfg, "results/figures"); figdir.mkdir(parents=True, exist_ok=True)
    for fn in (fig_phase_a, fig_phase_b, fig_abm):
        print(f"[ok] {fn(cfg, figdir)}")


if __name__ == "__main__":
    main()
