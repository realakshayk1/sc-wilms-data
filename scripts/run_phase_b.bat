@echo off
setlocal
cd /d "%~dp0.."
set WILMS_DEMO=0
echo === Phase B: Visium H and E histology ML ===
python phase2_histology_ml\01_extract_tiles.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\02_segment_nuclei.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\03_nucleus_features.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\04_train_classifier.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\05_spot_fractions.py --force %*
if errorlevel 1 exit /b 1
echo === Phase B validation ===
python phase2_histology_ml\08_loto_validation.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\09_negative_controls.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\10_snrna_spatial_concordance.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\11_marker_deconv_validation.py --force %*
if errorlevel 1 exit /b 1
python phase2_histology_ml\07_figures.py
echo === Phase B complete (PhysiCell excluded) ===
