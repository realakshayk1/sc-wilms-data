@echo off
REM Windows helper — R is often not on PATH after winget install.
set "RSCRIPT=C:\Program Files (x86)\R\R-4.6.1\bin\x64\Rscript.exe"
if not exist "%RSCRIPT%" set "RSCRIPT=C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
if not exist "%RSCRIPT%" (
  echo Rscript not found. Install: winget install RProject.R
  exit /b 1
)
cd /d "%~dp0.."
"%RSCRIPT%" %*
