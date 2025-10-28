#!/bin/bash
set -e
echo "🚀 Installing Julia and GRAPE capsule..."

# Installer Julia (Linux x86_64)
curl -fsSL https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.2-linux-x86_64.tar.gz -o julia.tar.gz
tar -xzf julia.tar.gz
export PATH="$PWD/julia-1.11.2/bin:$PATH"

# Cloner ton dépôt
git clone https://github.com/<ton_user>/GRAPE.git
cd GRAPE/capsule_light

# Initialiser l’environnement Julia et exécuter
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/example_ParkerSolarProbe.jl
