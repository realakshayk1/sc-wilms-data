@echo off
setlocal
cd /d "%~dp0.."
set WILMS_DEMO=0
echo === Phase B: Visium H&E histology ML ===
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
python phase2_histology_ml\06_map_to_physicell.py %*
python phase2_histology_ml\07_figures.py
echo === Phase B complete ===
