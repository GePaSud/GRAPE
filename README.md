# 🪐 GRAPE — General Relativity Accelerometer-based Propagation Environment

**Authors:** JP. Barriot, J. O’Leary, J. Ya, JM. Mari
**License:** MIT  
**Version:** 1.0 — October 2025  
**Journal:** *Software X* (submitted)

---

## 🚀 Overview

**GRAPE** (General Relativity Accelerometer-based Propagation Environment)  
is a Julia-based framework for simulating spacecraft trajectories in a **fully relativistic formulation**.

It integrates the motion of a spacecraft within arbitrary spacetime metrics (Schwarzschild, Kerr, Newtonian approximations, etc.), including **non-gravitational forces** and **accelerometer-based models**.

This repository contains:
- The **core engine** written in Julia (`src/GRAPE_core.jl`)
- A simplified **example** (Parker Solar Probe)
- Several **reproducible capsules** (Software X standard)

---

## 🧩 Capsule Types

GRAPE provides four reproducible environments:

| Capsule | Description | Reproducibility level |
|----------|--------------|------------------------|
| 🟩 **Light Capsule** | For users with Julia already installed | Code-only |
| 🟨 **Docker-Light Capsule** | Automatic build in official Julia Docker | Environment-fixed |
| 🟦 **Full Docker Capsule** | Includes menu and execution options | Interactive |
| 🟥 **Installer Capsule** | Installs Julia and runs GRAPE automatically (Windows/Linux) | Fully automated |

---

## 📁 Repository Structure
```text
GRAPE/
│
├── src/                            # Core physical and mathematical routines
│   └── GRAPE_core.jl
│
├── examples/                       # Demonstration cases
│   ├── example_ParkerSolarProbe.jl
│   └── example_simplified_no_tetrad.jl
│
├── capsule_julia_light/            # Light capsule (Julia already installed)
│   └── readme                      # Instructions for quick local run
│
├── capsule_docker_light/           # Docker light capsule (Julia official image)
│   ├── Dockerfile
│   └── readme                      # Minimal instructions
│
├── capsule_docker_full/            # Full interactive Docker capsule
│   ├── Dockerfile
│   ├── start_menu.sh               # Interactive shell menu (dialog-based)
│   └── readme                      # Usage instructions
│
├── capsule_julia_install/          # Auto-installers (Windows & Linux)
│   ├── install_grape.ps1           # PowerShell installer for Windows
│   ├── install_grape.sh            # Bash installer for Linux/macOS
│   └── readme                      # Usage instructions
│
├── LICENSE                         # MIT license
└── README.md                       # Main documentation (this file)
