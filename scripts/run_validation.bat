@echo off
setlocal
cd /d "%~dp0.."
echo === Phase A validation ===
call scripts\rscript.bat phase1_mechanotypes\09_distributional_validation.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\10_label_sensitivity.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\11_pseudobulk_validation.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\06_waddR_decompose.R || echo [warn] waddR skipped
call scripts\rscript.bat phase1_mechanotypes\08_figures.R || exit /b 1
echo === Phase B validation (requires Phase B pipeline outputs) ===
python phase2_histology_ml\08_loto_validation.py --force || exit /b 1
python phase2_histology_ml\09_negative_controls.py --force || exit /b 1
python phase2_histology_ml\10_snrna_spatial_concordance.py --force || exit /b 1
python phase2_histology_ml\11_marker_deconv_validation.py --force || exit /b 1
python phase2_histology_ml\05_spot_fractions.py --force || exit /b 1
python phase2_histology_ml\07_figures.py
echo === Validation complete ===
