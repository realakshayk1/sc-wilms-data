#!/usr/bin/env python3
"""WS1–WS4 visualization: render the coupling / bounding results into one dashboard + panels.

Reads the CSVs produced by coupling_core.R / 19_bifurcation_transfer.R / 20_bulk_coupling.R /
26_tissue_architecture.py and 05_uq.py, and draws:
  (1) tumor-level coupling network (FDR edges bold, red=neg/blue=pos)
  (2) tumor correlation heatmap
  (3) bifurcation — proliferation (bimodal) vs a unimodal lever
  (4) virtual-cohort coupling preservation (proliferation vs derived crowding)
  (5) sc-vs-bulk coupling recovery (what cell resolution buys)
  (6) tissue architecture across tumors (nodular vs diffuse)

Output: results/figures/phase_c_{dashboard, network, heatmap, forest, ...}.png (flat, phase-prefixed)
Usage: python phase1_mechanotypes/21_coupling_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CD = ROOT / "results" / "couplings"
OUT = ROOT / "results" / "figures"
PFX = "phase_c_"                                     # flat, phase-prefixed (matches phase_a_/phase_b_)
NEG, POS, GREY = "#c0392b", "#2471a3", "#95a5a6"     # colorblind-safe red / blue / grey
SHORT = {"proliferation": "prolif", "tp53_target": "p53", "wnt_canonical": "wnt",
         "blastemal_nephrogenic": "blastemal", "igf": "igf", "emt_axis": "EMT",
         "crowding_sensitivity": "crowding", "hypoxia_tolerance": "hypoxia"}


def _abbr(s):
    return SHORT.get(s, s)


def panel_network(ax):
    ed = pd.read_csv(CD / "network_tumorB_partial.csv")
    nodes = sorted(set(ed.a) | set(ed.b))
    ang = {n: 2 * np.pi * i / len(nodes) + np.pi / 2 for i, n in enumerate(nodes)}
    pos = {n: (np.cos(a), np.sin(a)) for n, a in ang.items()}
    for _, r in ed.iterrows():
        sig = r.bh_fdr < 0.10
        if not sig and abs(r.partial_r) < 0.25:
            continue
        x0, y0 = pos[r.a]; x1, y1 = pos[r.b]
        ax.plot([x0, x1], [y0, y1], color=(NEG if r.partial_r < 0 else POS),
                lw=(0.6 + 5 * abs(r.partial_r)) if sig else 0.7,
                alpha=(0.9 if sig else 0.18), zorder=1, solid_capstyle="round")
    for n, (x, y) in pos.items():
        ax.scatter([x], [y], s=760, color="white", edgecolor="#2c3e50", lw=1.5, zorder=2)
        ax.text(x, y, _abbr(n), ha="center", va="center", fontsize=8.5, zorder=3, weight="bold")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.4, 1.4); ax.axis("off")
    ax.set_title("Tumor-level coupling network\n(bold = BH-FDR<0.10)", fontsize=10, weight="bold")
    ax.legend(handles=[Line2D([0], [0], color=POS, lw=3, label="positive"),
                       Line2D([0], [0], color=NEG, lw=3, label="negative")],
              loc="lower center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.04))


def panel_heatmap(ax):
    m = pd.read_csv(CD / "network_tumorB_marginal_matrix.csv", index_col=0)
    m = m.rename(index=_abbr, columns=_abbr)
    im = ax.imshow(m.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(m))); ax.set_xticklabels(m.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index, fontsize=7)
    for i in range(len(m)):
        for j in range(len(m)):
            v = m.to_numpy()[i, j]
            if i != j:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 0.5 else "#2c3e50")
    ax.set_title("Program correlations (tumor level)", fontsize=10, weight="bold")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=7)


def panel_bifurcation(ax):
    ts = pd.read_csv(CD / "tumor_scores.csv")
    bif = pd.read_csv(CD / "bifurcation.csv").set_index("program")
    for prog, col, lab in [("proliferation", POS, "proliferation (bimodal)"),
                           ("blastemal_nephrogenic", GREY, "blastemal (unimodal)")]:
        ax.hist(ts[prog], bins=14, alpha=0.55, color=col, label=lab, edgecolor="white", lw=0.5)
    sp = bif.loc["proliferation", "split"]
    ax.axvline(sp, color=NEG, ls="--", lw=1.5, label=f"prolif split z={sp:.2f}")
    ax.set_xlabel("program score (z)", fontsize=8); ax.set_ylabel("tumors", fontsize=8)
    ax.set_title("Bifurcation: only proliferation splits the cohort", fontsize=10, weight="bold")
    ax.legend(fontsize=7, frameon=False); ax.tick_params(labelsize=7)


def panel_virtual(ax):
    d = pd.read_csv(ROOT / "results" / "abm" / "virtual_cohort" / "draws.csv")
    ax.scatter(d.proliferation, d.crowding_z, s=10, color=POS, alpha=0.45, edgecolor="none")
    r = np.corrcoef(d.proliferation, d.crowding_z)[0, 1]
    b = np.polyfit(d.proliferation, d.crowding_z, 1)
    xs = np.linspace(d.proliferation.min(), d.proliferation.max(), 20)
    ax.plot(xs, np.polyval(b, xs), color=NEG, lw=2)
    ax.set_xlabel("drawn proliferation (z)", fontsize=8)
    ax.set_ylabel("derived crowding (z)", fontsize=8)
    ax.set_title(f"Virtual cohort preserves coupling\n(prolif–crowding r={r:.2f}, n={len(d)})",
                 fontsize=10, weight="bold")
    ax.tick_params(labelsize=7)


def panel_recovery(ax):
    cmp = pd.read_csv(CD / "bulk" / "network_concordance.csv")
    cmap = {"recovered_by_bulk": POS, "NEEDS_cell_resolution": NEG,
            "bulk_only(composition?)": "#e67e22", "neither": GREY}
    for v, c in cmap.items():
        s = cmp[cmp.verdict == v]
        if len(s):
            ax.scatter(s.partial_r_sc, s.partial_r_bulk, s=70, color=c, label=v,
                       edgecolor="#2c3e50", lw=0.6, zorder=3)
    lim = 0.9
    ax.plot([-lim, lim], [-lim, lim], color=GREY, ls=":", lw=1, zorder=1)
    ax.axhline(0, color="#ddd", lw=0.8); ax.axvline(0, color="#ddd", lw=0.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("sc (cell-scoped) partial r", fontsize=8)
    ax.set_ylabel("bulk partial r", fontsize=8)
    ax.set_title("What cell resolution buys (coupling recovery)", fontsize=10, weight="bold")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left"); ax.tick_params(labelsize=7)


def panel_forest(ax):
    tr = pd.read_csv(CD / "transfer.csv")
    yt, ylab, y = [], [], 0.0
    for axis, g in tr.groupby("axis"):
        for _, r in g.sort_values("beta").iterrows():
            sig = r.boot_p < 0.05
            c = (NEG if r.beta < 0 else POS) if sig else GREY
            ax.plot([r.ci_lo, r.ci_hi], [y, y], color=c, lw=2.4,
                    alpha=0.95 if sig else 0.5, solid_capstyle="round", zorder=2)
            ax.scatter([r.beta], [y], color=c, s=42, zorder=3, edgecolor="white", lw=0.7)
            yt.append(y); ylab.append(f"{_abbr(axis)} ← {_abbr(r.lever)}{'*' if sig else ''}")
            y += 1
        y += 0.8
    ax.axvline(0, color="#888", ls="--", lw=1, zorder=1)
    ax.set_yticks(yt); ax.set_yticklabels(ylab, fontsize=7)
    ax.set_xlabel("standardized transfer β (95% CI)   * boot p<0.05", fontsize=8)
    ax.set_title("Intrinsic → extrinsic transfer rules", fontsize=10, weight="bold")
    ax.tick_params(labelsize=7); ax.invert_yaxis()


def panel_architecture(ax):
    a = pd.read_csv(ROOT / "results" / "spatial" / "tissue_architecture_per_sample.csv")
    a = a.sort_values("nodularity_index", ascending=True)
    y = np.arange(len(a))
    ax.barh(y, a.nodularity_index, color=POS, alpha=0.85)
    ax.set_yticks([]); ax.set_xlabel("nodularity index", fontsize=8)
    ax.set_ylabel(f"{len(a)} tumors (sorted)", fontsize=8)
    ax.set_title("Tissue architecture: nodular ↔ diffuse", fontsize=10, weight="bold")
    ax.tick_params(labelsize=7)


CLASS_COLOR = {"expressive": NEG, "intermediate": GREY, "sensitive": POS}
CLASS_ORDER = ["sensitive", "intermediate", "expressive"]


def sensitive_expressive_figure():
    tum = pd.read_csv(CD / "sensitive_expressive_tumor.csv")
    stats = pd.read_csv(CD / "sensitive_expressive_stats.csv")
    cell = CD / "sensitive_expressive_cellsub.csv"
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # A — intrinsic plane: proliferation vs p53, labeled by class
    ax = axes[0]
    for c in CLASS_ORDER:
        s = tum[tum["class"] == c]
        ax.scatter(s.proliferation, s.tp53_target, s=70, color=CLASS_COLOR[c],
                   edgecolor="#2c3e50", lw=0.6, label=f"{c} (n={len(s)})", zorder=3)
    ax.axvline(tum.proliferation.median(), color="#bbb", ls="--", lw=1)
    ax.axhline(tum.tp53_target.median(), color="#bbb", ls="--", lw=1)
    ax.text(0.98, 0.02, "expressive\n(hi prolif / lo p53)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=NEG, weight="bold")
    ax.text(0.02, 0.98, "sensitive\n(lo prolif / hi p53)", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color=POS, weight="bold")
    ax.set_xlabel("proliferation (z)", fontsize=9); ax.set_ylabel("p53 activity (z)", fontsize=9)
    ax.set_title("Intrinsic classification of 40 tumors", fontsize=10, weight="bold")
    ax.legend(fontsize=8, frameon=False, loc="lower left"); ax.tick_params(labelsize=8)

    # B — extrinsic axes by class (tumor level): the conditional distributions
    ax = axes[1]
    axes_names = ["crowding_sensitivity", "emt_axis", "hypoxia_tolerance"]
    for gi, axis in enumerate(axes_names):
        for ci, c in enumerate(CLASS_ORDER):
            vals = tum.loc[tum["class"] == c, axis].to_numpy()
            x = gi * 3 + ci
            ax.scatter(np.full(len(vals), x) + np.random.default_rng(0).normal(0, 0.06, len(vals)),
                       vals, color=CLASS_COLOR[c], s=28, alpha=0.8, edgecolor="none", zorder=3)
            ax.plot([x - 0.28, x + 0.28], [np.median(vals)] * 2, color="#2c3e50", lw=2, zorder=4)
        p = stats[(stats.level == "tumor") & (stats.axis == axis)].mwu_p.iloc[0]
        ax.text(gi * 3 + 1, 0.98, f"p={p:.1e}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8,
                weight="bold" if p < 0.05 else "normal")
    ax.axhline(0, color="#ddd", lw=0.8)
    ax.set_xticks([1, 4, 7]); ax.set_xticklabels(["crowding", "EMT", "hypoxia"], fontsize=9)
    ax.set_ylabel("extrinsic axis (z)", fontsize=9)
    ax.set_title("Extrinsic axes conditional on intrinsic class (tumor)", fontsize=10, weight="bold")
    ax.tick_params(labelsize=8)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=CLASS_COLOR[c],
                              markersize=8, label=c) for c in CLASS_ORDER], fontsize=8, frameon=False)

    # C — cell-level crowding by class: attenuated (couplings are between-tumor)
    ax = axes[2]
    if cell.exists():
        cs = pd.read_csv(cell)
        for ci, c in enumerate(CLASS_ORDER):
            v = cs.loc[cs["class"] == c, "crowding_sensitivity"].to_numpy()
            parts = ax.violinplot([v], positions=[ci], showmedians=True, widths=0.8)
            for b in parts["bodies"]:
                b.set_facecolor(CLASS_COLOR[c]); b.set_alpha(0.6)
        d = stats[(stats.level == "cell") & (stats.axis == "crowding_sensitivity")].cliffs_delta.iloc[0]
        ax.set_xticks(range(3)); ax.set_xticklabels(CLASS_ORDER, fontsize=8)
        ax.set_ylabel("cell crowding (z)", fontsize=9)
        ax.set_title(f"Same split at CELL level: near-flat\n(Cliff's δ={d:.2f} — couplings are between-tumor)",
                     fontsize=10, weight="bold")
        ax.tick_params(labelsize=8)
    fig.suptitle("Sensitive vs expressive: intrinsic state → extrinsic microenvironment",
                 fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / f"{PFX}sensitive_expressive.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    panels = [("network", panel_network), ("heatmap", panel_heatmap),
              ("forest", panel_forest), ("bifurcation", panel_bifurcation),
              ("virtual_cohort", panel_virtual), ("coupling_recovery", panel_recovery),
              ("architecture", panel_architecture)]
    fig, axes = plt.subplots(2, 4, figsize=(21, 9.5))
    flat = axes.ravel()
    for (name, fn), ax in zip(panels, flat):
        try:
            fn(ax)
        except Exception as e:  # pragma: no cover
            ax.text(0.5, 0.5, f"{name}:\n{e}", ha="center", va="center", fontsize=8, wrap=True)
            ax.axis("off")
    flat[len(panels)].axis("off")   # 8th cell (legend space) blank
    fig.suptitle("Wilms coupled-lever & resolution-bounding results (SCPCP000006, n=40)",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / f"{PFX}dashboard.png", dpi=150, bbox_inches="tight")
    # also save each panel standalone
    for name, fn in panels:
        f, a = plt.subplots(figsize=(6, 5))
        try:
            fn(a)
            f.tight_layout(); f.savefig(OUT / f"{PFX}{name}.png", dpi=150, bbox_inches="tight")
        except Exception:  # pragma: no cover
            pass
        plt.close(f)
    try:
        sensitive_expressive_figure()
        print("[ok] sensitive_expressive.png")
    except Exception as e:  # pragma: no cover
        print(f"[warn] sensitive_expressive figure skipped: {e}")
    print(f"[ok] figures -> {OUT}/{PFX}dashboard.png (+ panels + {PFX}sensitive_expressive)")


if __name__ == "__main__":
    main()
