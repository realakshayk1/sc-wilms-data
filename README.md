# sc-wilms-data

**Wilms tumor distributional mechanotypes + histology-informed spatial composition**

Computational pipeline connecting the Radhakrishnan lab's **Wasserstein mechanotyping framework** (Phase A) to **Visium H&E morphology ML** (Phase B) for the public ScPCA cohort [**SCPCP000006**](https://scpca.alexslemonade.org/projects/SCPCP000006) — paired snRNA-seq, Visium spatial transcriptomics, H&E images, and bulk RNA-seq from favorable and anaplastic Wilms tumors.

> **Scientific goal:** Identify how whole *distributions* of interpretable molecular programs differ across Wilms compartments and histology, then validate whether H&E-derived spatial cell-state composition agrees with transcriptomics — a prerequisite for morphology-informed agent-based models (PhysiCell).

---

## Table of contents

1. [Background](#background)
2. [Pipeline overview](#pipeline-overview)
3. [Data & cohort](#data--cohort)
4. [Phase A — Mechanotypes (snRNA-seq)](#phase-a--mechanotypes-snrna-seq)
5. [Phase B — Histology ML (Visium H&E)](#phase-b--histology-ml-visium-he)
6. [Results summary](#results-summary)
7. [Methods tried & what they taught us](#methods-tried--what-they-taught-us)
8. [ABM initial conditions](#abm-initial-conditions-physicell)
9. [Figure gallery](#figure-gallery)
10. [Quick start](#quick-start)
11. [Repository layout](#repository-layout)
12. [Configuration](#configuration)
13. [Limitations & next steps](#limitations--next-steps)
14. [References & citation](#references--citation)

---

## Background

Wilms tumor (nephroblastoma) is morphologically organized into **blastemal**, **epithelial**, and **stromal** compartments, with clinically dominant **favorable vs anaplastic** histology. Phase A asks *where* the favorable/anaplastic and relapse signal lives by interrogating the snRNA-seq with several lenses — the lab's **distributional mechanotype** (whole-distribution shifts via **Wasserstein-1** + **consensus clustering**), compartment **composition**, moderated **differential expression**, and **pathway enrichment** — and reporting the signal wherever it actually appears. For Wilms the answer is **compositional + pathway-level** (a relapse-associated proliferation program), not a within-compartment distributional shift.

Separately, spatial ABM models often assume uniform or deconvolution-only initial conditions. Phase B tests what paired **H&E** can recover: it robustly reads **anaplasia** (the prognostically decisive nuclear-atypia phenotype) but **not** continuous compartment composition — so H&E contributes the tumor's growth *regime*, while compartment fractions for the ABM come from the transcriptomic deconvolution.

---

## Pipeline overview

```mermaid
flowchart TB
  subgraph data [ScPCA SCPCP000006]
    SN[snRNA-seq 40 samples]
    VIS[Visium + H&E 41 samples]
  end

  subgraph phaseA [Phase A — Mechanotypes R]
    QC[QC + cell-state mapping]
    SC[1-D program scores]
    W1[Wasserstein-1 per feature]
    CC[ConsensusClusterPlus]
    SW[Mechanotype switches]
    QC --> SC --> W1 --> CC --> SW
  end

  subgraph phaseB [Phase B — Histology ML Python]
    TIL[H&E spot tiles]
    EMB[Phikon-v2 embeddings]
    SEG[StarDist nuclei + morphology]
    CLS[Anaplasia classifier - MIL]
    TIL --> EMB --> CLS
    TIL --> SEG --> CLS
  end

  subgraph abm [ABM initial conditions]
    PC[Per-tumor PhysiCell parameters]
  end

  SN --> QC
  VIS --> TIL
  CC --> PC
  CLS --> PC
```

**Reproducibility:** All stochastic steps use seed `42` (logged). Intermediates live in `data/processed/`; headline outputs in `results/`.

---

## Data & cohort

| Modality | Count | Use |
|----------|-------|-----|
| snRNA-seq (nucleus) | 40 samples | Phase A mechanotypes |
| Visium spots + H&E | 41 samples | Phase B morphology ML |
| Bulk RNA-seq | 45 samples | Optional validation (not wired) |

**Histology:** 23 favorable / 22 anaplastic (`subdiagnosis` in ScPCA metadata).

**Access:**
- **Metadata (no token):** `python scripts/fetch_scpca_metadata.py`
- **API download:** `ScPCAr` R package ([docs](https://alexslemonade.github.io/ScPCAr/))
- **Manual download (recommended on Windows):** Portal zips → `scripts/ingest_manual_downloads.ps1`

Raw data are **never committed**; provenance logged in `data/raw/scpca_access_log.txt`.

---

## Phase A — Mechanotypes (snRNA-seq)

### Methodology

| Step | Script | Method |
|------|--------|--------|
| Ingest | `scripts/ingest_manual_scpca.R` | Load merged SCE; join `subdiagnosis` histology |
| QC + labels | `02_qc_normalize.R` | ≥200 genes/cell; assign compartments from **fetal-kidney signatures** (CM/UB/PV/fibroblast, `config/cell_signatures.yaml`) on tumor cells — gene-wise z-scored module scores + top-vs-runner-up margin gate |
| Scores | `03_compute_scores.R` | **Fixed** gene programs (`config/features.yaml`): log1p(CPM<sub>pos</sub>) − log1p(CPM<sub>neg</sub>) via `gene_symbol` |
| Items | `04_wasserstein_matrix.R` | Groups = (compartment × histology); **≥25 cells** rule |
| Distance | `04_wasserstein_matrix.R` | **1-D Wasserstein-1 only** on score distributions (`transport` package) |
| Clustering | `05_consensus_cluster.R` | ConsensusClusterPlus PAM; k via low **PAC** + high **Calinski–Harabasz** |
| Switches | `07_mechanotype_switches.R` | Flag compartment if cluster assignment differs favorable vs anaplastic |
| Composition | `12_composition_analysis.R` | Per-sample compartment fractions; CLR-Wilcoxon, patient-level, BH-FDR |
| Moderated DE | `14_moderated_de.R` | Pseudobulk **edgeR-QLF** + **limma-voom**, histology & relapse axes |
| Pathway GSEA | `15_hallmark_gsea.R` | **fgsea** preranked on the limma-voom moderated *t*, 50 MSigDB Hallmark sets |
| Prognostics | `16_prognostic_association.R` | **Firth** logistic + Fisher of composition/proliferation vs relapse; bootstrap/profile CI |
| Figures | `08_figures.R`, `18_result_figures.py` | W1/switch heatmaps, score violins, PAC/CHI, GSEA + DE summary |

Methods log: `results/mechanotypes/phase_a_methods.yaml`

### Key design choices

- **No feature fishing:** the program list is predefined before clustering (blastemal, epithelial, stromal, proliferation, WT1, Wnt/β-catenin); sensitivity to it is reported rather than tuned.
- **Compartment labels from fetal-kidney biology:** reference/`cellassign` annotations call WT tumor cells "hemangioblast/trophoblast/Unknown" and cannot separate the three compartments, so compartments are assigned from **fetal-kidney developmental signatures** on tumor cells (z-scored module scores + margin gate).
- **Patient is the unit of inference:** histology/relapse contrasts are tested across the ~40 samples, never across cells (cell-level testing is pseudoreplication and inflates significance).
- **1-D Wasserstein only:** multivariate Wasserstein on gene matrices underperforms on scRNA-seq (lab benchmark); W1 runs only on the predefined 1-D scores.
- **Match the instrument to the biology:** the same data is interrogated by distributional, compositional, single-gene DE, and pathway tests — and the signal is reported wherever it actually lives (composition + pathway here), with method-robust negatives reported as negatives.

### What the Wasserstein mechanotype measures

Most single-cell comparisons collapse a cell group to a **mean** or a **proportion**. The
mechanotype approach instead compares the **whole distribution** of a program's activity. For two
1-D score distributions, the **Wasserstein-1 (earth-mover's) distance** is the area between their
cumulative distribution functions — the minimum "work" to reshape one into the other — so it
responds to differences in **location, spread, and shape** that a difference-in-means misses
(two groups can share a mean yet differ in how bimodal or heavy-tailed a program is).

Per program, the pairwise W1 matrix over (compartment × histology) groups is **consensus-clustered**
into distributional classes — the "mechanotypes" — and a compartment is said to **switch
mechanotype** if its class differs between favorable and anaplastic. W1 is applied **only to the
predefined 1-D scores**: multivariate Wasserstein is statistically unstable on high-dimensional
scRNA-seq (lab benchmark), so the high-dim gene matrix is never passed to it.

**Outcome for Wilms — a method-robust negative.** With valid (fetal-kidney) labels and
patient-level statistics, **no** compartment shows a within-compartment distributional shift
(0/18 at BH-FDR, both axes). The machinery is sound; the phenomenon simply isn't there — the
favorable/anaplastic and relapse signal lives in **compartment composition** and a **relapse
proliferation program**, neither of which a *within-compartment distribution* comparison can see.
That negative is exactly what motivated the composition / moderated-DE / pathway analyses below.

### Phase A results

Compartments are assigned from **fetal-kidney developmental signatures** (cap mesenchyme,
ureteric bud, primitive vesicle, fibroblast; `config/cell_signatures.yaml`) on tumor cells.
All inference is **patient-level** — labels are contrasted across the ~40 samples, never across
cells (cell-level testing is pseudoreplication) — with BH-FDR within each analysis.

Wilms histology separates along two molecular dimensions:

**1 · Compartment composition** (the histology axis). The *relative abundance* of compartments
differs by histology: epithelial fraction ↑ in anaplastic (0.59 vs 0.44, BH-p=0.004), the
PV/mature-epithelial subgroup ↑ in anaplastic (BH-p=0.005), stromal ↑ in favorable (BH-p=0.038).
Within-compartment program *distributions*, by contrast, do not differ (0/18 at BH-FDR<0.05 on
both axes): the histology signal lives in composition, not in shifted distributions of program
activity.

**2 · A proliferative transcriptional program** (the relapse axis), resolved at three levels:

![Phase A — Hallmark GSEA and moderated DE](results/figures/phase_a_gsea_de.png)

| Analysis | Script | Result |
|----------|--------|--------|
| **Moderated DE** | `14_moderated_de.R` | edgeR-QLF: **130 genes FDR<0.05** for histology (NOTCH2, PODXL, PTPRO, DACT3), 39 for relapse |
| **Hallmark GSEA** | `15_hallmark_gsea.R` | **166 significant pathway-contrasts**; on the relapse axis **E2F_TARGETS (q=9e-29)**, **G2M_CHECKPOINT**, **MYC_TARGETS** are up, replicated in the epithelial (q=4e-28) and stromal (q=6e-30) compartments |
| **Prognostics** | `16_prognostic_association.R` | A pseudobulk **proliferation score predicts relapse** (Firth OR≈4/SD, p=0.013; Fisher OR 7.6, p=0.017) — nominal (not BH-FDR-significant; n=10 relapse) |

The relapse axis — **higher cell-cycle/E2F/G2M/MYC activity** — is the canonical aggressive-tumor
program and matches the literature ([Yang 2025](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1539897/full); TP53/anaplasia). Gene-level DE,
pathway GSEA, and patient-level prognostics converge on it independently. Overall survival is not
modelable in this cohort (`vital_status` has only 5 deaths).

Details: `composition_analysis.csv`, `moderated_de.csv`, `hallmark_gsea.csv`, `prognostic_association.csv`

---

## Phase B — Histology ML (Visium H&E)

### Methodology

| Step | Script | Method |
|------|--------|--------|
| Tiles | `01_extract_tiles.py` | Visium hires H&E patches centered on tissue spots; **Macenko** stain norm (ref `SCPCL000438`) |
| Programs | (in 01) | Same Phase A gene scores on spot RNA → dominant state + softmax fractions |
| Segment | `02_segment_nuclei.py` | **StarDist** `2D_versatile_he` (learned H&E nuclei model) |
| Features | `03_nucleus_features.py` | Area, eccentricity, solidity, texture, H-intensity, neighbor density |
| Embeddings | `15_phase_b_mil.py` | **Phikon-v2** (ViT-L, 1024-d) tile embeddings, 200 spots/tumor |
| Classifier | `14`/`15_phase_b_mil.py` | Histology (anaplasia) from embeddings — mean-pool + **attention-MIL**, leave-one-tumor-out |
| Morphology | `16_stardist_morphology.py` | Per-tumor nuclear-atypia features (giant-nucleus fraction, pleomorphism) → RF; embedding ensemble |
| Composition | `12_spot_composition_regression.py`, `13_fm_embedding_regression.py` | Cross-modal LOTO regression: predict per-spot compartment fractions from morphology / FM embeddings (negative, with shuffled + random controls) |
| Stats | `phase_b_stats.py` | **DeLong** AUC CIs, label-permutation p, paired DeLong |
| ABM | `06`/`17_positives_to_abm.py` | Map composition + proliferation + anaplasia → per-tumor PhysiCell parameters |
| Figures | `07_figures.py`, `18_result_figures.py` | Deconv scatter, AUC forest, ABM parameter panels |

Methods log: `results/classifier/phase_b_methods.json` · Config: `config/phase_b.yaml`

### Phase B results

Full cohort: **41 tumors, ~260k Visium spots, H&E at hires resolution.**

**Reading anaplasia from H&E.** Unfavorable histology is *defined* by nuclear atypia (giant,
hyperchromatic, pleomorphic nuclei; [Vujanić 2024](https://onlinelibrary.wiley.com/doi/full/10.1002/pbc.31000)) — the signal H&E carries
natively. Tumor-level classification of anaplastic vs favorable, held out across tumors
(leave-one-tumor-out), with DeLong 95% CIs and label-permutation p:

![Phase B — histology AUC forest](results/figures/phase_b_histology_auc.png)

| Model | AUC | 95% CI (DeLong) | perm p |
|-------|-----|-----------------|--------|
| **Phikon-v2 embeddings, attention-MIL** | **0.748** | [0.57, 0.87] | 0.003 |
| Phikon-v2 embeddings, mean-pool | 0.733 | [0.55, 0.86] | 0.006 |
| StarDist nuclear morphology | 0.687 | [0.50, 0.83] | 0.021 |
| Ensemble (morphology + embedding) | 0.719 | [0.53, 0.85] | 0.009 |

All models predict anaplasia significantly above chance. The signal **saturates near AUC
~0.73–0.75**: attention-MIL is statistically indistinguishable from flat mean-pooling (paired
DeLong p=0.83), and the morphology+embedding ensemble does not exceed the embedding alone
(p=0.57) — both read the same nuclear atypia, bounded by Visium-hires tile resolution (median ~14
segmentable nuclei/tumor). Learned StarDist segmentation is essential: classical watershed
morphology reaches only AUC 0.39.

**H&E does not read continuous composition.** The complementary cross-modal task — predicting
per-spot transcriptomic compartment *fractions* from morphology — is a clean negative:
leave-one-tumor-out Pearson *r* ≈ 0 for both hand-crafted features and FM embeddings (≈
shuffled/random controls). H&E therefore sets the tumor's **growth regime** (anaplastic ⇒
aggressive); fine compartment composition for ABM initial conditions comes from the transcriptomic
deconvolution.

Details: `phase_b_mil_phikon-v2.json`, `stardist_morphology.json`, `fm_embedding_regression_phikon.json`

---

## Results summary

| Finding | Result |
|---------|--------|
| **Composition shifts by histology** | Epithelial ↑ anaplastic, stromal ↑ favorable — BH-FDR<0.05 (CLR-Wilcoxon) |
| **Proliferation program marks relapse** | E2F/G2M/MYC GSEA q≈1e-29; edgeR-QLF 130 DE genes FDR<0.05; proliferation score → relapse (Firth p=0.013) |
| **Within-compartment distributions** | No shift (0/18 BH-FDR) — histology lives in composition, not distribution |
| **H&E predicts anaplasia** | Tumor-level AUC **0.748** (attention-MIL), perm p=0.003; StarDist morphology 0.687 (p=0.021) |
| **H&E predicts continuous composition** | No — LOTO *r*≈0 (FM + hand-crafted); H&E sets growth regime, not fine composition |
| **ABM initial conditions** | Per-tumor PhysiCell parameters from composition + proliferation + anaplasia |
| Reproducible repo | Pinned env, numbered scripts, config-driven paths, DeLong/permutation/Firth stats, unit tests |

---

## Methods tried & what they taught us

Several approaches were tried and ruled out. Each negative is **method-robust** (held out,
FDR/permutation-corrected) and pointed to the analysis that did work — they are findings, not
dead ends.

![Methods that didn't work](results/figures/methods_negatives.png)

*Left:* every within-compartment Wasserstein test sits far below the BH-FDR=0.05 line (0/36).
*Middle:* the Welch t-test finds ~0 single-gene hits where edgeR-QLF moderation finds up to 130.
*Right:* H&E→compartment-composition regression gives held-out *r*≈0, indistinguishable from
shuffled/random controls. (The watershed-segmentation, scale/MIL, and ensemble negatives are in
the [Phase B AUC forest](#phase-b--histology-ml-visium-he).)

| Approach | What it tested | Outcome | What it taught us |
|----------|----------------|---------|-------------------|
| **Wasserstein-1 distributional mechanotype** | within-compartment shifts in the *distribution* of program activity | ✗ 0/18 BH-FDR (both axes) | the signal is **compositional + pathway-level**, not a within-compartment distribution shift (explained above) |
| **Reference / `cellassign` cell labels** | compartment identity | ✗ invalid for WT | reference annotations call tumor cells "hemangioblast/trophoblast/Unknown" → assign from **fetal-kidney signatures** instead |
| **Cell-level permutation testing** | contrast significance | ✗ anti-conservative (every p≈0.001) | cells within a tumor aren't independent (pseudoreplication) → permute labels **across patients** |
| **Welch t-test pseudobulk DE** | single-gene differential expression | ✗ ~0 genes at FDR<0.05 | underpowered at n≈20/group → empirical-Bayes moderation (**edgeR-QLF / limma-voom**) recovers 130 |
| **Watershed nuclei segmentation** | nuclear morphology → histology | ✗ AUC 0.39 (worse than chance) | classical segmentation is too weak on hires tiles → **StarDist** lifts it to 0.687 |
| **H&E → continuous composition regression** (hand-crafted + Phikon) | read 3-compartment fractions from morphology | ✗ leave-one-tumor-out *r*≈0 | H&E carries the tumor's **regime** (anaplasia), not fine composition; resolution-bound |
| **Scale spots + ViT-L encoder + attention-MIL** | lift the histology AUC | ✗ 0.724→0.748, MIL vs mean-pool p=0.83 | ~0.73 is a **Visium-hires resolution ceiling**, not a modeling gap |
| **Morphology + embedding ensemble** | add signal orthogonal to the embedding | ✗ no gain (paired p=0.57) | both modalities read the **same** nuclear atypia |

What survived these is in [Results summary](#results-summary): compartment composition, the
relapse proliferation program (DE + GSEA + prognostics), H&E anaplasia classification, and the
StarDist morphology fix.

---

## ABM initial conditions (PhysiCell)

The Phase A/B findings translate directly into per-tumor agent-based-model inputs
([`17_positives_to_abm.py`](phase2_histology_ml/17_positives_to_abm.py) →
[`results/abm/positives_to_physicell.yaml`](results/abm/positives_to_physicell.yaml)):

![ABM parameters](results/figures/abm_parameters.png)

| Finding | PhysiCell parameter |
|---------|---------------------|
| Compartment composition | **initial cell-type fractions** (right panel) |
| Proliferation score | **proliferation_rate** multiplier, bounded `1 + 0.6·z` |
| p53-target activity | **apoptosis_rate** multiplier |
| H&E anaplasia probability | **high-grade regime** (extra proliferation, reduced adhesion) |

The mapping encodes the measured biology rather than being fit to outcome: as a check, the
proliferation multiplier averages **1.40 in relapse vs 0.98 in non-relapse** (left panel). This
gives PhysiCell spatially-resolved, biologically-grounded starting conditions per tumor rather than
a uniform configuration. Running the simulation itself requires the PhysiCell binary (cluster).

---

## Figure gallery

Regenerate all figures:

```powershell
scripts\run_figures.bat
```

| Figure | Description |
|--------|-------------|
| [`phase_a_w1_heatmaps.png`](results/figures/phase_a_w1_heatmaps.png) | Pairwise W1 distances per feature program |
| [`phase_a_mechanotype_switch_heatmap.png`](results/figures/phase_a_mechanotype_switch_heatmap.png) | Switches favorable ↔ anaplastic by compartment |
| [`phase_a_score_distributions.png`](results/figures/phase_a_score_distributions.png) | WT1, blastemal, proliferation score violins |
| [`phase_a_consensus_metrics.png`](results/figures/phase_a_consensus_metrics.png) | PAC & CHI vs k for consensus clustering |
| [`phase_b_deconv_validation.png`](results/figures/phase_b_deconv_validation.png) | H&E vs RNA spot fractions (3 compartments) |
| [`phase_b_dominant_state_confusion.png`](results/figures/phase_b_dominant_state_confusion.png) | Dominant compartment agreement matrix |
| [`phase_b_fractions_by_histology.png`](results/figures/phase_b_fractions_by_histology.png) | Composition by favorable vs anaplastic |
| [`phase_b_segmentation_mosaic.png`](results/figures/phase_b_segmentation_mosaic.png) | Segmentation QC on sample tiles |
| [`phase_b_classifier_summary.png`](results/figures/phase_b_classifier_summary.png) | Accuracy metrics + correlation bar chart |
| [`mechanotype_switches.png`](results/figures/mechanotype_switches.png) | Bar chart of switches (from script 07) |
| [`phase_a_gsea_de.png`](results/figures/phase_a_gsea_de.png) | Hallmark GSEA (relapse axis) + moderated-DE FDR gene counts |
| [`phase_b_histology_auc.png`](results/figures/phase_b_histology_auc.png) | Histology AUC forest with DeLong 95% CIs (watershed→StarDist→Phikon→MIL→ensemble) |
| [`abm_parameters.png`](results/figures/abm_parameters.png) | ABM proliferation multiplier by relapse + per-tumor initial fractions |
| [`methods_negatives.png`](results/figures/methods_negatives.png) | Method-robust negatives: distributional mechanotype, Welch DE, H&E→composition |

Phase A/B figures regenerate with `08_figures.R` + `python phase2_histology_ml/18_result_figures.py`.
Segmentation overlays: `data/processed/nuclei/overlays/`

---

## Quick start

### 1. Environment

```powershell
# Python deps (Phase B)
pip install scanpy scikit-learn scikit-image opencv-python-headless pyarrow pyyaml matplotlib seaborn scipy

# R 4.x + packages (Phase A)
winget install RProject.R
scripts\rscript.bat scripts\install_r_packages.R
scripts\rscript.bat scripts\scpca_auth.R   # optional for API download
```

Or: `conda env create -f environment.yml && conda activate sc-wilms-data`

### 2. Metadata (no download)

```powershell
python scripts/fetch_scpca_metadata.py
```

### 3. Manual data ingest (recommended)

Place Portal zips in Downloads, then:

```powershell
powershell -File scripts/ingest_manual_downloads.ps1
scripts\rscript.bat scripts\ingest_manual_scpca.R
```

### 4. Run pipelines

```powershell
# Phase A: QC → mechanotypes → figures
scripts\run_phase_a.bat

# Phase B: requires spaceranger extract; sets WILMS_DEMO=0
scripts\run_phase_b.bat

# Figures only
scripts\run_figures.bat
```

### 5. Tests

```powershell
pytest -q
scripts\rscript.bat -e "testthat::test_dir('tests', filter = 'phase1')"
```

---

## Repository layout

```
sc-wilms-data/
├── config/                  # paths.yaml, features.yaml, phase_b.yaml, physicell.yaml
├── phase1_mechanotypes/     # R: 00–16 numbered scripts (mechanotype, composition, DE, GSEA, prognostics)
├── phase2_histology_ml/     # Python: 00–18 numbered scripts (tiles, embeddings, MIL, StarDist, ABM, figures)
├── scripts/                 # ingest, run_phase_*.bat, fetch metadata
├── data/raw/                # gitignored ScPCA downloads
├── data/processed/          # gitignored intermediates (SCE, scores, tiles, nuclei, embeddings)
├── results/
│   ├── mechanotypes/        # composition, moderated DE, GSEA, prognostics CSVs; consensus RDS
│   ├── classifier/          # histology AUC + MIL/StarDist JSON, deconv/composition JSON
│   ├── figures/             # analysis figures (PNG)
│   └── abm/                 # per-tumor PhysiCell parameters (YAML + CSV)
├── tests/
├── PRD.md                   # requirements & acceptance criteria
├── AGENTS.md                # agent/human coding rules
└── .learnings/LEARNINGS.md  # accumulated gotchas
```

---

## Configuration

| File | Purpose |
|------|---------|
| `config/features.yaml` | **Fixed** Phase A gene programs, consensus params, seed |
| `config/paths.yaml` | All relative paths (no hard-coded absolutes) |
| `config/phase_b.yaml` | Visium tile size, library limits, segmentation, classifier |
| `config/physicell.yaml` | ABM domain and cell-type parameter mapping |

**Important:** Set `WILMS_DEMO=0` (or use `run_phase_b.bat`) for real Visium processing. `WILMS_DEMO=1` generates synthetic tiles for CI only.

---

## Limitations & next steps

1. **Compartment assignment** uses fetal-kidney signatures on tumor cells; validate against spatial/IHC ground truth where available.
2. **Phase A coverage:** only ~30% of nuclei map confidently to a compartment after the margin gate.
3. **Segmentation:** StarDist `2D_versatile_he` needs a Windows directory *junction* (not a symlink — avoids the admin requirement) to load.
4. **waddR decomposition** (`06_waddR_decompose.R`): optional location/shape/size interpretation.
5. **PhysiCell:** initial-condition mapping is produced (`positives_to_physicell.yaml`); the simulation itself needs the PhysiCell binary on a cluster.

### Externally-gated extensions

Two ceilings are set by the available data, not by method, and need inputs beyond the local
ScPCA cohort to lift:

- **Phase B resolution.** Only Visium-**hires** tiles (~96 px/spot) are available, not the original
  whole-slide images — capping tumor-level anaplasia AUC near 0.73 and StarDist near ~14 nuclei/tumor.
  Lifting it needs the raw WSIs + a gated pathology foundation model (**UNI2 / Virchow2 /
  Prov-GigaPath**) or **XMAG** (5×-native), gated on an HF access token.
- **Time-to-event survival.** Local metadata is **binary** (`relapse_status`) with too few deaths
  (`vital_status`, n=5) for Cox / Kaplan-Meier. Proper recurrence-free-survival validation of the
  proliferation signature needs **TARGET-WT** (GDC) or GSE31403/GSE10320, plus a **Scissor**
  reproduction of the relapse-cell analysis.

### Status against the project goal

- **Which compartments shift transcriptional behavior?** Resolved: the shift is **compositional**
  (epithelial ↑ anaplastic, FDR<0.05) plus a **proliferation program on the relapse axis**
  (E2F/G2M/MYC q≈1e-29; 130 moderated-DE genes). Within-compartment program *distributions* do not
  shift — a clean negative that localizes the signal to composition.
- **Does H&E track composition well enough to seed the ABM?** H&E robustly reads **anaplasia**
  (AUC ~0.73, held out) — the prognostically decisive feature — but **not** continuous
  3-compartment composition (cross-tumor *r*≈0). H&E therefore sets the tumor's **growth regime**;
  fine composition comes from transcriptomic deconvolution. The ~0.73 is resolution-bound (above).
- **(3) — yes, as a mapping.** Per-tumor PhysiCell parameters are generated and pass a directional
  sanity check. What remains is **running PhysiCell itself** (the binary, on a cluster) and ideally
  cross-cohort survival validation (Tier-3).

**Bottom line:** the analysis half of the project is *complete and honestly characterized* for this
cohort — every locally-answerable question has an answer with effect sizes, CIs, and stated nulls.
The two open ends are **external** (higher-resolution histology; time-to-event survival) and the
**downstream PhysiCell simulation**, none of which are blocked by missing analysis code.

---

## References & citation

- ScPCA Portal & `SCPCP000006`: [Alex's Lemonade ScPCA](https://scpca.alexslemonade.org/) · preprint [10.1101/2024.04.19.590243](https://doi.org/10.1101/2024.04.19.590243)
- ScPCAr R package: [GitHub](https://github.com/AlexsLemonade/ScPCAr)
- Wasserstein mechanotyping framework: Radhakrishnan lab pan-cancer mechanobiology work (W1 + ConsensusClusterPlus + n≥25 rule)
- waddR: 2-Wasserstein decomposition for scRNA-seq
- Visium: 10x Genomics spatial; Macenko stain normalization; StarDist H&E nuclei (`2D_versatile_he`)

**Lab context:** Multiscale intrinsic–extrinsic coupling (PhysiCell ABM) — see `PRD.md` appendix.

---

## ScPCAr API cheat sheet

| Task | Function | Auth? |
|------|----------|-------|
| List projects | `scpca_projects()` | No |
| Sample table | `get_project_samples("SCPCP000006")` | No |
| Agree to terms | `get_auth(email, agree=TRUE)` | — |
| Download merged SCE | `download_project(..., format="sce", merged=TRUE)` | Yes |
| Download Visium | `download_project(..., format="spatial")` | Yes |

**Deprecated:** `download_sample()` — use `create_dataset()` → `download_dataset(await_processing=TRUE)`.

See also [ScPCA download guide](https://scpca.readthedocs.io/en/stable/download_files.html).
