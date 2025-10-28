#!/usr/bin/env julia
# ============================================================
#  GRAPE Example — Parker Solar Probe trajectory simulation
#  Compatible with module GRAPE_P
#  Author: Jean-Pierre Barriot, October 2025
# ============================================================

# --- Load GRAPE core code
include("../src/GRAPE_core.jl")
using .GRAPE_P
using Dates, Printf, LinearAlgebra, PythonPlot

# --- Output directory
mkpath("output")

# --- Constants (SI units)
const c_light = BigFloat(2.99792458e8)          # speed of light (m/s)
const mu      = BigFloat(1.32712440018e20)      # Sun GM (m^3/s^2)
const rs      = BigFloat(2*mu/(c_light^2))      # Schwarzschild radius
const as      = BigFloat(0.0)                   # rotation (a=0 for Schwarzschild)
const AU      = BigFloat(1.495978707e11)        # astronomical unit (m)
const solar_radiation_pressure = BigFloat(4.56e-6)  # at 1 AU (N/m^2)
const mass_PSP = BigFloat(655.0)                # spacecraft mass (kg)
const shield_area = BigFloat(4.0)               # m^2
const shield_reflection_coef = BigFloat(1.8)

# --- Integration parameters
const evpar_step = BigFloat(10.0)               # step in seconds
const N_integration_steps = 2000                # number of steps
const time_option = 1                           # proper time integration
const order_parder_metric = 4                   # metric derivative order
const itmax_rki = 8                             # fixed-point iterations
const epsmin_CSN = BigFloat(1e-6)
const epsabs_CSN = BigFloat(1e-6)
const eps_tetrad = BigFloat(1e-10)

# --- Initial conditions (Parker Solar Probe at perihelion)
r0 = BigFloat(9.86 * 6.957e8)                   # ≈ 9.86 R☉
θ0 = BigFloat(pi/2)
ϕ0 = BigFloat(0)
v0_tan = sqrt(mu / r0) * BigFloat(0.6)          # 60% circular
# The GRAPE state uses alternating (x,v) for 4D coordinates: (t, r, φ, θ)
# Initialize time, radius, etc. in polar coordinates
y0 = fill(BigFloat(0), 26)
y0[1] = BigFloat(0)                             # coordinate time
y0[3] = r0
y0[5] = ϕ0
y0[7] = θ0
# initial 4-velocity (approximate, non-relativistic)
y0[2] = BigFloat(1.0)                           # dt/dτ (normalized later)
y0[4] = BigFloat(0.0)
y0[6] = v0_tan
y0[8] = BigFloat(0.0)
y0[9] = BigFloat(0.0)                           # proper time
y0[10] = BigFloat(0.0)                          # integration time
# tetrad initialization (identity)
for i in 1:4
    for j in 1:4
        y0[4*i+j+6] = i == j ? BigFloat(1.0) : BigFloat(0.0)
    end
end

# --- Prepare initial conditions in Cholesky frame
ychol0 = GRAPE_P.ychol_from_ychart(y0, GRAPE_P.Kerr_Metric_Polar)

# --- Run integration
@printf("Launching Parker Solar Probe test case...\n")
timetag = Dates.format(now(), "yyyymmdd_HHMMSS")
(evparf, yfinal, dtaudevpar_internal) = GRAPE_P.main(
    BigFloat(0.0),
    ychol0,
    GRAPE_P.Kerr_Metric_Polar
)
@printf("Integration complete. Final evpar = %.2e s\n", Float64(evparf))

# --- Read results from output file
data = readdlm("SCephemeris_$(timetag).txt")
idx = reshape(data[:, 1], (26, :))
vals = reshape(data[:, 2], (26, :))
r_vals = vals[3, :]
ϕ_vals = vals[5, :]

x = r_vals .* cos.(ϕ_vals)
y = r_vals .* sin.(ϕ_vals)

# --- Plot trajectory
figure(figsize=(6,6))
plot(x ./ (6.957e8), y ./ (6.957e8), lw=1.5)
xlabel("x / R☉")
ylabel("y / R☉")
title("Parker Solar Probe orbit (GRAPE relativistic test)")
axis("equal")
grid(true)
savefig("output/PSP_orbit.png", dpi=200)
@printf("Orbit plot saved to output/PSP_orbit.png\n")

println("✅ Test completed successfully.")
