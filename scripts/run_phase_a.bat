@echo off
cd /d "%~dp0.."
call scripts\rscript.bat phase1_mechanotypes\02_qc_normalize.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\03_compute_scores.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\04_wasserstein_matrix.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\05_consensus_cluster.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\07_mechanotype_switches.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\06_waddR_decompose.R || echo [warn] waddR skipped
call scripts\rscript.bat phase1_mechanotypes\09_distributional_validation.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\10_label_sensitivity.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\11_pseudobulk_validation.R || exit /b 1
call scripts\rscript.bat phase1_mechanotypes\08_figures.R || exit /b 1
echo [ok] Phase A complete. See results/mechanotypes/ and results/figures/
