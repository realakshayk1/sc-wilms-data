@echo off
REM Run full demo pipeline (Windows). Set WILMS_DEMO=1 for synthetic data.
set WILMS_DEMO=1
cd /d "%~dp0.."

echo === Phase A (R demo) ===
Rscript phase1_mechanotypes\01_download.R --demo
Rscript phase1_mechanotypes\02_qc_normalize.R
Rscript phase1_mechanotypes\03_compute_scores.R
Rscript phase1_mechanotypes\04_wasserstein_matrix.R
Rscript phase1_mechanotypes\05_consensus_cluster.R
Rscript phase1_mechanotypes\07_mechanotype_switches.R

echo === Phase B (Python demo) ===
python phase2_histology_ml\01_extract_tiles.py --demo
python phase2_histology_ml\02_segment_nuclei.py --demo
python phase2_histology_ml\03_nucleus_features.py
python phase2_histology_ml\04_train_classifier.py
python phase2_histology_ml\05_spot_fractions.py
python phase2_histology_ml\06_map_to_physicell.py

echo === Done ===
