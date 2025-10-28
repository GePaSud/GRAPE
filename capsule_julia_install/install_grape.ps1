Write-Host "🚀 Installation complète de GRAPE (Windows)..." -ForegroundColor Cyan

# 1. Vérifier si Julia est installée
if (-not (Get-Command julia -ErrorAction SilentlyContinue)) {
    Write-Host "🔧 Julia non trouvée. Téléchargement..." -ForegroundColor Yellow
    $url = "https://julialang-s3.julialang.org/bin/winnt/x64/1.11/julia-1.11.2-win64.exe"
    $file = "$env:TEMP\julia-installer.exe"
    Invoke-WebRequest $url -OutFile $file
    Start-Process -FilePath $file -ArgumentList "/SILENT", "/DIR=$env:USERPROFILE\Julia-1.11.2" -Wait
    $env:Path += ";$env:USERPROFILE\Julia-1.11.2\bin"
}

# 2. Cloner ton dépôt GitHub
if (-not (Test-Path "$env:USERPROFILE\GRAPE")) {
    Write-Host "📦 Clonage du dépôt GRAPE..." -ForegroundColor Cyan
    git clone https://github.com/<ton_user>/GRAPE.git "$env:USERPROFILE\GRAPE"
}

# 3. Aller dans la capsule light et exécuter
cd "$env:USERPROFILE\GRAPE\capsule_light"
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. examples/example_ParkerSolarProbe.jl
