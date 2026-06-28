@echo off
setlocal
cd /d "%~dp0.."
echo === Regenerating all figures ===
call scripts\rscript.bat phase1_mechanotypes\08_figures.R || exit /b 1
set WILMS_DEMO=0
py -3 phase2_histology_ml\07_figures.py 2>nul || python phase2_histology_ml\07_figures.py || exit /b 1
echo === Figures written to results/figures/ ===
