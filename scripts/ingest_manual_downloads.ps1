# Ingest manually downloaded ScPCA zips into data/raw/scpca_downloads/
param(
    [string]$DownloadsDir = "$env:USERPROFILE\Downloads",
    [string]$AccessDate = "2026-06-28",
    [switch]$Wait,
    [int]$PollSeconds = 30
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path $PSScriptRoot -Parent
$RawDir = Join-Path $RepoRoot "data\raw"
$ZipArchive = Join-Path $RawDir "scpca_downloads\zips"
$ExtractRoot = Join-Path $RawDir "scpca_downloads"
New-Item -ItemType Directory -Force -Path $ZipArchive, $ExtractRoot | Out-Null

$Expected = @(
    "SCPCP000006_single-cell-experiment_$AccessDate.zip",
    "SCPCP000006_spaceranger_$AccessDate.zip"
)

function Get-ZipIfReady($name) {
    $path = Join-Path $DownloadsDir $name
    if (-not (Test-Path $path)) { return $null }
    $size1 = (Get-Item $path).Length
    Start-Sleep -Seconds 3
    $size2 = (Get-Item $path).Length
    if ($size1 -ne $size2) { return $null }  # still downloading
    return $path
}

function Wait-ForZips {
    Write-Host "Waiting for downloads in $DownloadsDir ..."
    while ($true) {
        $ready = @()
        foreach ($n in $Expected) {
            $p = Get-ZipIfReady $n
            if ($p) { $ready += $p } else { Write-Host "  pending: $n" }
        }
        if ($ready.Count -eq $Expected.Count) { return $ready }
        Start-Sleep -Seconds $PollSeconds
    }
}

if ($Wait) {
    $sources = Wait-ForZips
} else {
    $sources = @()
    foreach ($n in $Expected) {
        $p = Join-Path $DownloadsDir $n
        if (Test-Path $p) {
            $sources += $p
        } else {
            Write-Warning "Not found: $p"
        }
    }
}

if (-not $sources -or $sources.Count -eq 0) {
    Write-Host @"

No zip files ready. When browser finishes, run:

  scripts\ingest_manual_downloads.ps1 -Wait

Or if files are elsewhere:

  scripts\ingest_manual_downloads.ps1 -DownloadsDir 'C:\path\to\folder'
"@
    exit 1
}

foreach ($src in $sources) {
    $dest = Join-Path $ZipArchive (Split-Path $src -Leaf)
    if (-not (Test-Path $dest) -or (Get-Item $src).Length -ne (Get-Item $dest).Length) {
        Write-Host "[move] $src -> $dest"
        Copy-Item $src $dest -Force
    }
    $kind = if ($dest -match "single-cell") { "single-cell-experiment" } else { "spaceranger" }
    $out = Join-Path $ExtractRoot $kind
    $marker = Join-Path $out ".extracted"
    New-Item -ItemType Directory -Force -Path $out | Out-Null
    if (-not (Test-Path $marker)) {
        Write-Host "[unzip] $dest -> $out (this may take several minutes)..."
        Expand-Archive -Path $dest -DestinationPath $out -Force
        Set-Content $marker (Get-Date -Format o)
    } else {
        Write-Host "[skip] already extracted: $out"
    }
}

$log = Join-Path $RawDir "scpca_access_log.txt"
Add-Content $log "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') | manual Portal download access_date=$AccessDate"

Write-Host "[ok] Ingest complete. Running R ingest..."
& (Join-Path $PSScriptRoot "rscript.bat") (Join-Path $PSScriptRoot "ingest_manual_scpca.R")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[ok] Starting Phase A pipeline..."
& (Join-Path $PSScriptRoot "run_phase_a.bat")
