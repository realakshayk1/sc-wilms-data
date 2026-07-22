#!/usr/bin/env python3
"""Second figure set for the coupling/bounding work:
  survival.png            — TARGET-WT OS: KM (blastemal, p53) + HR forest of all 9 levers
  virtual_cohort_ranges.png — the PhysiCell parameter ranges the sweep feeds (the deliverable)
  tissue_maps.png         — real Visium compartment maps: nodular vs diffuse tumor
  lever_heatmap.png       — 40 tumors x 8 levers, clustered, annotated by favorable/anaplastic
Output: results/figures/couplings/. Usage: python phase1_mechanotypes/28_more_figures.py
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
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "phase3_abm"))
from abm_utils import load_config, load_spot_coords_um, resolve_path  # noqa: E402

CD = ROOT / "results" / "couplings"
OUT = ROOT / "results" / "figures" / "couplings"
NEG, POS, GREY = "#c0392b", "#2471a3", "#95a5a6"
COMP = {"blastemal": "#8e44ad", "epithelial": "#16a085", "stromal": "#e67e22", "necrotic": "#7f8c8d"}
SHORT = {"proliferation": "prolif", "tp53_target": "p53", "wnt_canonical": "wnt",
         "blastemal_nephrogenic": "blastemal", "igf": "igf", "emt_axis": "EMT",
         "crowding_sensitivity": "crowding", "hypoxia_tolerance": "hypoxia"}


def _lev_label(s):
    return s.replace("lever_", "").replace("_hi", " (hi)").replace("_mut", " (mut)").replace(
        "blastemal_nephrogenic", "blastemal").replace("tp53_target", "p53").replace("wnt_canonical", "wnt")


def km_curve(time, event):
    time = np.asarray(time, float); event = np.asarray(event, int)
    xs, ys, s = [0.0], [1.0], 1.0
    for ut in np.unique(time):
        n = (time >= ut).sum(); d = ((time == ut) & (event == 1)).sum()
        if n > 0 and d > 0:
            s *= 1 - d / n
        xs.append(ut); ys.append(s)
    return np.array(xs), np.array(ys)


def fig_survival():
    sd = pd.read_csv(CD / "target_wt_surv_df.csv")
    res = pd.read_csv(CD / "survival_levers.csv")
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    for ax, (col, title) in zip(axes[:2],
                                [("lever_blastemal_nephrogenic_hi", "Blastemal program"),
                                 ("lever_tp53_target_hi", "p53-target activity")]):
        for val, color, lab in [(1, NEG, "high"), (0, POS, "low")]:
            g = sd[sd[col] == val]
            if len(g) < 3:
                continue
            x, y = km_curve(g.time, g.event)
            ax.step(x, y, where="post", color=color, lw=2.2,
                    label=f"{lab} (n={len(g)}, deaths={int(g.event.sum())})")
        row = res[res.lever == col]
        ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("OS probability")
        ax.set_title(f"{title}\nHR={row.hr.iloc[0]:.2f}  p={row.p_univariate.iloc[0]:.2f}",
                     fontsize=10.5, weight="bold")
        ax.legend(fontsize=8.5, frameon=False, loc="lower left")
    ax = axes[2]
    r = res.dropna(subset=["hr"]).sort_values("hr").reset_index(drop=True)
    for i, row in r.iterrows():
        c = NEG if row.hr > 1 else POS
        sig = row.p_univariate < 0.05
        ax.plot([row.hr_lo, row.hr_hi], [i, i], color=c, lw=2.2, alpha=0.9 if sig else 0.6,
                solid_capstyle="round")
        ax.scatter([row.hr], [i], color=c, s=48, zorder=3, edgecolor="white", lw=0.6)
    ax.axvline(1, color="#888", ls="--", lw=1)
    ax.set_yticks(range(len(r))); ax.set_yticklabels([_lev_label(l) for l in r.lever], fontsize=8)
    ax.set_xscale("log"); ax.set_xlabel("hazard ratio (OS, log scale) — 95% CI")
    ax.set_title("OS hazard ratios\n(0/9 survive BH-FDR<0.10; trends only)", fontsize=10.5, weight="bold")
    fig.suptitle("TARGET-WT overall survival (125 RNA / 38 MAF cases): honest trends, no FDR hits",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "survival.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_virtual_ranges():
    d = pd.read_csv(ROOT / "results" / "abm" / "virtual_cohort" / "draws.csv")
    params = [("proliferation_mult", "proliferation rate ×"), ("apoptosis_mult", "apoptosis rate ×"),
              ("adhesion_mult", "adhesion ×"), ("motility_mult", "migration speed ×"),
              ("igf_uptake_mult", "IGF2 uptake ×"), ("pressure_half_max", "pressure→cycle half-max"),
              ("necrosis_half_max", "O₂→necrosis half-max")]
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for (col, lab), ax in zip(params, axes.ravel()):
        v = d[col].to_numpy()
        ax.hist(v, bins=30, color=POS, alpha=0.8, edgecolor="white", lw=0.4)
        lo, med, hi = np.percentile(v, [5, 50, 95])
        ax.axvspan(lo, hi, color=NEG, alpha=0.10)
        ax.axvline(med, color=NEG, lw=2)
        ax.set_title(f"{lab}\nmedian {med:.2f} · p5–p95 [{lo:.2f}, {hi:.2f}]", fontsize=9.5, weight="bold")
        ax.set_yticks([]); ax.tick_params(labelsize=8)
    axes.ravel()[-1].axis("off")
    axes.ravel()[-1].text(0.5, 0.5, "These are the bounded, coupled\nranges the PhysiCell sensitivity\n"
                          "analysis samples from\n(256 correlated synthetic tumors).",
                          ha="center", va="center", fontsize=11, style="italic", color="#444")
    fig.suptitle("Virtual-cohort PhysiCell parameter ranges (the seeding deliverable)",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT / "virtual_cohort_ranges.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_tissue_maps():
    cfg = load_config()
    sig = pd.read_parquet(resolve_path(cfg, "data/processed/spot_signatures.parquet"))
    sig = sig[sig["in_tissue"] == 1]
    spot_um = float(cfg["phase_c"].get("spot_diameter_um", 55.0))
    root = resolve_path(cfg, str(cfg["paths"]["phase_b"]["spatial_root"]))
    lib_dirs = {p.name.replace("_spatial", ""): p for p in root.glob("SCPCS*/SCPCL*_spatial")}
    # both tumors are blastemal-rich (~40-50%) — the contrast is spatial ORGANIZATION, not composition
    picks = [("SCPCS000200", "Nodular blastemal (37% blast · autocorr 0.42)"),
             ("SCPCS000173", "Diffuse blastemal (50% blast · autocorr 0.07)")]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    for ax, (sid, title) in zip(axes, picks):
        g = sig[sig.sample_id == sid]
        libid = g.library_id.iloc[0]
        coords = load_spot_coords_um(lib_dirs[libid], spot_um)
        m = g.merge(coords[["x_um", "y_um"]], left_on="barcode", right_index=True, how="inner")
        for comp, c in COMP.items():
            s = m[m.dominant_state == comp]
            ax.scatter(s.x_um, s.y_um, s=9, color=c, label=f"{comp} ({len(s)})", edgecolor="none")
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{title}\n{sid} · {len(m)} spots", fontsize=11, weight="bold")
        ax.legend(fontsize=8, frameon=False, markerscale=2, loc="upper right")
    fig.suptitle("Real Visium compartment maps → ABM initial seeding geometry", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "tissue_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_lever_heatmap():
    ts = pd.read_csv(CD / "tumor_scores.csv")
    meta = pd.read_csv(CD / "tumor_meta.csv")
    levers = ["proliferation", "tp53_target", "wnt_canonical", "blastemal_nephrogenic",
              "igf", "emt_axis", "crowding_sensitivity", "hypoxia_tolerance"]
    M = ts.set_index("sample_id")[levers]
    ro = leaves_list(linkage(M.to_numpy(), method="ward"))          # cluster tumors
    co = leaves_list(linkage(M.to_numpy().T, method="ward"))        # cluster levers
    Mo = M.iloc[ro, co]
    hist = meta.set_index("sample_id")["subdiagnosis"].reindex(Mo.index).fillna("?")
    hcol = hist.str.lower().map(lambda s: POS if "favor" in s else (NEG if "anapl" in s else GREY))

    fig, (axc, ax) = plt.subplots(1, 2, figsize=(9, 11), gridspec_kw={"width_ratios": [0.04, 1]})
    axc.imshow(np.array([[0]] * len(Mo)), aspect="auto", cmap="Greys", vmin=0, vmax=1)
    for i, c in enumerate(hcol):
        axc.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, color=c))
    axc.set_xlim(-0.5, 0.5); axc.set_ylim(len(Mo) - 0.5, -0.5); axc.axis("off")
    axc.set_title("histology", fontsize=8)
    im = ax.imshow(Mo.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(co))); ax.set_xticklabels([SHORT[Mo.columns[i]] for i in range(len(co))],
                                                      rotation=45, ha="right", fontsize=9)
    ax.set_yticks([]); ax.set_ylabel(f"{len(Mo)} tumors (ward-clustered)", fontsize=9)
    ax.set_title("Per-tumor intrinsic-lever profiles", fontsize=11, weight="bold")
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="lever score (z)")
    fig.legend(handles=[Line2D([0], [0], marker="s", color="w", markerfacecolor=POS, markersize=9, label="favorable"),
                        Line2D([0], [0], marker="s", color="w", markerfacecolor=NEG, markersize=9, label="anaplastic")],
               loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(OUT / "lever_heatmap.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in [("survival", fig_survival), ("virtual_cohort_ranges", fig_virtual_ranges),
                     ("tissue_maps", fig_tissue_maps), ("lever_heatmap", fig_lever_heatmap)]:
        try:
            fn(); print(f"[ok] {name}.png")
        except Exception as e:
            import traceback; print(f"[FAIL] {name}: {e}"); traceback.print_exc()


if __name__ == "__main__":
    main()
