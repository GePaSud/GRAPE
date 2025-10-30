# ==================================================================================================
# 🌌 ROOT SECTION — GRAPE NUMERICAL CONSTANTS AND GAUSS 5-STAGE INTEGRATOR
# --------------------------------------------------------------------------------------------------
# These constants are the backbone of the GRAPE numerical environment.
# They are propagated throughout the entire codebase and used for all high-precision computations.
# Implementation adapted for Windows 11 (Julia BigFloat environment).
# ==================================================================================================

# --- Load GRAPE Core -------------------------------------------------------
include("../src/GRAPE_core.jl")
# ---------------------------------------------------------------------------
# ================================================================================================
# ⚙️ BASIC NUMERICAL CONSTANTS
# ------------------------------------------------------------------------------------------------
# Define all fundamental BigFloat constants used in integrator formulas.
# Each constant is defined explicitly as BigFloat to maintain numerical consistency.
# ================================================================================================
const BF0   = BigFloat(0)       ; const BF1   = BigFloat(1)       ; const BF6   = BigFloat(6)
const BF2   = BigFloat(2)       ; const BF8   = BigFloat(8)       ; const BF11  = BigFloat(11)
const BF13  = BigFloat(13)      ; const BF23  = BigFloat(23)
const BF32  = BigFloat(32)      ; const BF35  = BigFloat(35)      ; const BF59  = BigFloat(59)
const BF63  = BigFloat(63)      ; const BF64  = BigFloat(64)      ; const BF70  = BigFloat(70)
const BF225 = BigFloat(225)     ; const BF308 = BigFloat(308)     ; const BF322 = BigFloat(322)
const BF405 = BigFloat(405)     ; const BF452 = BigFloat(452)     ; const BF960 = BigFloat(960)
const BF1080 = BigFloat(1080)   ; const BF3240 = BigFloat(3240)   ; const BF3600 = BigFloat(3600)
const BF05 = BF1 / BF2
const BF3 = BigFloat(3) ; const BF4 = BigFloat(4) ; const BF5 = BigFloat(5) ; const BF12 = BigFloat(12) ; const BF20 = BigFloat(20)
const BF60 = BigFloat(60) ; const BF105 = BigFloat(105) ; const BF280 = BigFloat(280)	
const BF180 = BigFloat(180) ; const BF1000 = BigFloat(1000)

# ================================================================================================
# 🧮 GAUSS 5-STAGE IMPLICIT SYMPLECTIC INTEGRATOR (Butcher 1962)
# ------------------------------------------------------------------------------------------------
# The following constants define the coefficients of the 5-stage Gauss–Legendre collocation method,
# an implicit symplectic Runge–Kutta integrator. It preserves the Hamiltonian structure of
# relativistic trajectories with high stability and precision.
#
# Reference: J.C. Butcher, *Implicit Runge–Kutta Processes*, Math. Comp. 18 (1964), 50–64.
# ================================================================================================

# --- Coefficients definitions --------------------------------------------------
const om1 = (BF322 - BF13 * sqrt(BF70)) / BF3600
const op1 = (BF322 + BF13 * sqrt(BF70)) / BF3600

const om2 = BF05 * sqrt((BF35 + BF2 * sqrt(BF70)) / BF63)
const op2 = BF05 * sqrt((BF35 - BF2 * sqrt(BF70)) / BF63)

const om3 = om2 * (BF452 + BF59 * sqrt(BF70)) / BF3240
const op3 = op2 * (BF452 - BF59 * sqrt(BF70)) / BF3240

const om4 = om2 * (BF64 + BF11 * sqrt(BF70)) / BF1080
const op4 = op2 * (BF64 - BF11 * sqrt(BF70)) / BF1080

const om5 = BF8 * om2 * (BF23 - sqrt(BF70)) / BF405
const op5 = BF8 * op2 * (BF23 + sqrt(BF70)) / BF405

const om6 = om2 - BF2 * om3 - om5
const op6 = op2 - BF2 * op3 - op5

const om7 = om2 * (BF308 - BF23 * sqrt(BF70)) / BF960
const op7 = op2 * (BF308 + BF23 * sqrt(BF70)) / BF960

# ================================================================================================
# 📊 BUTCHER TABLEAUX (a_B, b_B, c_B)
# ------------------------------------------------------------------------------------------------
# - a_B : Coefficient matrix of the Runge–Kutta method (implicit coupling)
# - b_B : Weight vector for stage averaging
# - c_B : Abscissa vector (stage time positions)
#
# These are used by GRAPE’s core integrator to compute each step in proper time evolution.
# ================================================================================================
const a_B = [
    om1          op1 - om3 + op4   BF32 / BF225 - om5   op1 - om3 - op4   om1 - om6     ;  # 5×5 coefficient matrix
    om1 - om3 + om4   op1          BF32 / BF225 - op5   op1 - op6         om1 - om3 - om4 ;
    om1 + om7          op1 + op7    BF32 / BF225         op1 - op7         om1 - om7      ;
    om1 + op3 + op4    op1 + op6    BF32 / BF225 + op5   op1               om1 + op3 - op4;
    om1 + op6          op1 + op3 + op4   BF32 / BF225 + om5   op1 + op3 - op4   om1
]

const b_B = [BF2 * om1  BF2 * op1  BF64 / BF225  BF2 * op1  BF2 * om1]
# Line vector — weights for the stage combination (∑ bᵢ * fᵢ)

const c_B = [BF05 - om2 , BF05 - op2 , BF05 , BF05 + op2 , BF05 + om2]
# Column vector — abscissas for internal stages (time fractions)
# ================================================================================================
# ==================================================================================================
# 🌍 GLOBAL INITIALIZATION — GRAPE ENVIRONMENT SETUP
# --------------------------------------------------------------------------------------------------
# This section initializes the GRAPE runtime environment:
#  - imports dependencies (Dates, PythonPlot, etc.)
#  - creates a unique timestamped working directory for each run
#  - opens the main run log file
#  - records basic session information (time, version, working path)
# ==================================================================================================

# --- Load core Julia libraries ------------------------------------------------
using Dates          # Access to system clock for time tagging
using PythonPlot     # Plotting backend (can be replaced by another library if needed)
# using SPICE        # Optional: activate for DE430 planetary ephemerides (commented by default)
# ------------------------------------------------------------------------------

# ================================================================================================
# 🕓 TIMESTAMP AND OUTPUT DIRECTORY SETUP
# ------------------------------------------------------------------------------------------------
# Each GRAPE run generates a unique time tag (yyyymmddHHMMSS) used to:
#  - identify the run in log files and figures
#  - name the output directory and associated data files
# ================================================================================================
td = "$(now())"
timetag = td[1:4] * td[6:7] * td[9:10] * td[12:13] * td[15:16] * td[18:19]  # Format: YYYYMMDDHHMMSS

# Create and move into run-specific output directory
mkdir("EXGR_" * timetag)
cd("EXGR_" * timetag)

# ================================================================================================
# 📝 RUN LOG INITIALIZATION
# ------------------------------------------------------------------------------------------------
# The run log file stores execution details, diagnostic messages, and timing information.
# It allows later traceability between figures, computed data, and the original configuration.
# ================================================================================================
runlogfile = open("runlog_$(timetag).txt", "w")

prt3("(G)eneral (R)elativity (A)ccelerometer-based (P)ropagation (E)nvironment, version July 7, 2025")
prt3("timetag = ", timetag, " ==> all files generated in this run are labeled with this timetag")
prt3("and are output to the working directory below:")

# Determine the working directory for this run
working_directory = pwd()
prt3("working directory = ", working_directory)

# Record the system time at the start of the computation
t_local = now()
prt3("computer clock at start of run = ", t_local)
# ------------------------------------------------------------------------------------------------
# ================================================================================================
# ⚙️ SYSTEM SIZE DEFINITION
# ------------------------------------------------------------------------------------------------
# nsyst defines the number of coupled differential equations integrated by GRAPE:
#  - 8  : equations of motion (4 position + 4 velocity components)
#  - 2  : proper time and integration time
#  - 16 : tetrad equations (4×4 orthonormal frame)
# ================================================================================================
const nsyst = 26
# ------------------------------------------------------------------------------------------------
# ==================================================================================================
# 🧠 NUMERICAL PRECISION CONTROL — MPFR BigFloat CONFIGURATION
# --------------------------------------------------------------------------------------------------
# GRAPE operates natively with arbitrary-precision floating-point numbers (BigFloat),
# using the GNU-MPFR library built into Julia.
#
# This section:
#   - sets the working precision (in bits)
#   - prints diagnostic information to confirm accuracy
#   - compares Julia’s internal π value with a high-precision reference
#
# The goal is to ensure full control over numerical precision across all GRAPE computations.
# ==================================================================================================

# ================================================================================================
# ⚙️ Set precision for BigFloat calculations
# ------------------------------------------------------------------------------------------------
# Common values for IEEE754 emulation:
#   - 53 bits  →  double precision (64-bit)
#   - 113 bits →  quadruple precision (128-bit)
#   - 237 bits →  extended precision (~256-bit)
#
# This precision setting applies globally to all BigFloat operations.
# ================================================================================================
setprecision(BigFloat, 237)

# Retrieve and log the active precision level
precision_bits = precision(BigFloat)
prt3("BigFloats bits for MPFR library = ", precision_bits)

# ================================================================================================
# 🔍 Validate numerical accuracy with π
# ------------------------------------------------------------------------------------------------
# Compare the current MPFR-computed π with an external reference value.
# This test verifies that the precision setting and MPFR environment are correctly configured.
# ================================================================================================
prt3("control of accuracy of MPFR library with respect to external value of PI ==>")

# π value from the active MPFR context
pi_run = string(one(BigFloat) * π)
prt3(pi_run, length(pi_run) - 1, "digits for this run from MPFR library")

# High-precision external reference (for control)
pi_ref = "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342"
prt3(pi_ref, length(pi_ref) - 1, "digits of pi for control")
# -----------------------------------------------------------------------------------------------


# ================================================================================================
# ⚠️ Important note on BigFloat conversion
# ------------------------------------------------------------------------------------------------
# Avoid direct conversions like:
#     BigFloat(2.1)  or  2.1 * BF1
# These produce rounding artifacts because 2.1 is stored as a binary float before conversion.
#
# Use instead:
#     parse(BigFloat, "2.1")
# to ensure exact decimal representation.
#
# This warning is especially important for Julia versions < 1.9.
# ================================================================================================
# NEVER do this type of stuff in the code: BigFloat(2.1) or 2.1*BF1, except for π which is auto-adaptive.
# Example of bad behavior (for Julia < 1.9):
#     2.100000000000000088817841970012523233890533447265625000000000000000000000000000
# ================================================================================================

#--------------------------------------------------------------------------------------------------------------------------------------------------
# ==================================================================================================
# 🧮 INTEGRATION AND NUMERICAL CONTROL PARAMETERS
# --------------------------------------------------------------------------------------------------
# This section defines all the numerical constants that control the behavior and accuracy
# of GRAPE’s implicit symplectic Runge–Kutta integrator.
#
# These parameters influence:
#   - the convergence of the fixed-point solver,
#   - the precision of metric derivatives (Christoffel symbols),
#   - the tolerance for tetrad orthonormality tests,
#   - the choice and step of the integration time variable.
# ==================================================================================================

# ================================================================================================
# 🔁 IMPLICIT RUNGE–KUTTA ITERATION CONTROL
# ------------------------------------------------------------------------------------------------
# GRAPE uses a symplectic implicit Runge–Kutta method based on fixed-point iterations.
#   - The number of iterations (itmax_rki) controls convergence accuracy.
#   - A higher value ensures stability in highly elliptical or relativistic orbits.
# ================================================================================================
const itmax_rki = 10  # Recommended ≥ 6; 10 is needed for Parker Solar Probe
prt3("number of iterations for implicit integrator = ", itmax_rki)

# This Runge–Kutta scheme is *symplectic*:
# it preserves the norm of the 4-velocity, ensuring energy conservation in relativistic space-time dynamics.
# ================================================================================================

# ================================================================================================
# 🧭 METRIC DERIVATIVE ACCURACY
# ------------------------------------------------------------------------------------------------
# The variable order_parder_metric defines the order of precision used for the
# numerical computation of metric derivatives (Christoffel symbols).
# Valid choices are: 2, 4, 6, or 8 (typically 4 for good balance).
# ================================================================================================
order_parder_metric = 4
prt3("central difference order for Christoffel symbols = ", order_parder_metric)
# ================================================================================================

# ================================================================================================
# ⚙️ NUMERICAL EPSILONS FOR METRIC DERIVATIVES
# ------------------------------------------------------------------------------------------------
# epsabs_CSN and epsmin_CSN define the denominator epsilon in central difference
# approximations when computing Christoffel symbol derivatives.
#
#   eps = max(epsmin_CSN, epsabs_CSN * abs(coordinate))
#
# These values should be adjusted if the coordinate scales vary significantly
# (e.g., near black holes or in extreme relativistic regimes).
# ================================================================================================
const epsabs_CSN = parse(BigFloat, "1.e-12")  # absolute scaling factor
const epsmin_CSN = parse(BigFloat, "1.e-50")  # lower limit for epsilon
prt3("first constant for numerical derivation of Christoffel symbols = ", Float64(epsabs_CSN))
prt3("second constant for numerical derivation of Christoffel symbols = ", Float64(epsmin_CSN))
# ================================================================================================

# ================================================================================================
# 🧱 TETRAD VALIDATION TOLERANCE
# ------------------------------------------------------------------------------------------------
# eps_tetrad controls the threshold for checking the orthonormality of the
# Fermi–Walker transported tetrad (formula 1 in the GRAPE reference paper).
#
# Large deviations indicate a breakdown of numerical stability.
# ================================================================================================
const eps_tetrad = parse(BigFloat, "1.e-9")
prt3("threshold for the numerical test of the validity of Fermi-Walker tetrad = ", Float64(eps_tetrad))
# ================================================================================================

# ================================================================================================
# ⏱️ INTEGRATION TIME DEFINITION
# ------------------------------------------------------------------------------------------------
# The variable time_option selects which notion of time is used for integration:
#   1 → proper time τ
#   2 → user-defined function of proper time
#   3 → coordinate time t (default for most runs)
#   4 → user-defined function of coordinate time
# The mapping between these times is handled in evpar2tau_external().
# ================================================================================================
const time_option = 3
prt3("integration time option = ", time_option)
prt3("1 ==> proper time, 2 ==> user defined time wrt proper time, 3 ==> coordinate time, 4 ==> user defined time wrt coordinate time")
# ================================================================================================

# ================================================================================================
# ⌛ INTEGRATION STEP AND DURATION
# ------------------------------------------------------------------------------------------------
# evpar_step: time step (in seconds) for the integration variable.
# N_integration_steps: total number of iterations to perform.
# The total simulated duration = N_integration_steps × evpar_step.
# ================================================================================================
const evpar_step = parse(BigFloat, "120.00")  # step in seconds
prt3("integration time step = ", Float64(evpar_step), " seconds")

const evpar0 = BF0  # starting value for integration parameter
prt3("integration time at start of integration = ", Float64(evpar0), " seconds")

N_integration_steps = 1_000  # typical run length (100_000 for long simulations)
prt3("number of loops for integration process = ", N_integration_steps)
prt3("duration period of probe orbit in this run = ", Float64(N_integration_steps * evpar_step), " seconds of integration time")
# ================================================================================================


# ================================================================================================
# 🌠 PHYSICAL CONSTANTS (auxiliary)
# ------------------------------------------------------------------------------------------------
# Speed of light, fixed by the IAU definition, expressed in km/s.
# Note: conversion to meters is handled elsewhere as needed.
# ================================================================================================
c_light_nature = parse(BigFloat, "299792.458")  # km·s⁻¹
# prt3("c_light_in_nature = ", Float64(c_light_nature), " km s⁻¹")
# ================================================================================================
# ==================================================================================================
# 🌌 PHYSICAL CONSTANTS AND SPACETIME METRIC DEFINITION
# --------------------------------------------------------------------------------------------------
# This section sets the fundamental constants and selects the metric model (Kerr, Schwarzschild, etc.)
# that defines the geometry of spacetime in which the spacecraft propagates.
#
# Constants are defined in *consistent relativistic units* (kilometers, seconds, kilograms),
# and all calculations are performed with arbitrary precision (BigFloat).
# ==================================================================================================


# ================================================================================================
# ⚡ SPEED OF LIGHT
# ------------------------------------------------------------------------------------------------
# The value of c is defined by the IAU standard (in km/s).
# It can be rescaled (e.g., divided by 10) for testing numerical sensitivity.
# ================================================================================================
const c_light = c_light_nature  # natural constant value; reduced for test cases if needed
prt3("c light = ", Float64(c_light), " km s^-1")

# Reciprocal of c² — used in metric and potential expansions.
# Setting eps_c2 = 0 reproduces the flat Minkowski metric.
const eps_c2 = BF1 / c_light^2
# ================================================================================================


# ================================================================================================
# 🌀 SPACETIME METRIC SELECTION
# ------------------------------------------------------------------------------------------------
# The active metric defines how curvature affects motion:
#   - Kerr_Metric_Polar → rotating (Kerr) spacetime in Boyer–Lindquist coordinates
#   - Schwarzschild_Metric_Polar → non-rotating (Schwarzschild) spacetime
#
# The routine Give_Function_Name() returns the symbolic name of the metric function,
# used for logging and runtime reference.
# ================================================================================================
const Kerr_Metric_Polar_name = Give_Function_Name(Kerr_Metric_Polar)  # required for runtime linkage
# Alternative example (for static Schwarzschild spacetime):
# Given_Metric = Give_Function_Name(Schwarzschild_Metric_Polar)
# ================================================================================================


# ================================================================================================
# 🌠 GRAVITATIONAL CONSTANT AND SOLAR PARAMETERS
# ------------------------------------------------------------------------------------------------
# G — Newtonian gravitational constant (km³·kg⁻¹·s⁻²)
# μ_s — gravitational parameter of the Sun (G × M_sun) in km³·s⁻²
# AU — astronomical unit, standard average Sun–Earth distance (km)
# ================================================================================================
const G = parse(BigFloat, "6.6743e-20")  # gravitational constant [km³/kg/s²]
# prt3("constant of gravitation G = ", Float64(G), " km³ kg⁻¹ s⁻²")

AU = parse(BigFloat, "149597870.700")   # 1 astronomical unit [km]
# prt3("one astronomical unit = ", Float64(AU), " km")

const mu_s = parse(BigFloat, "0.132712440018e12")  # solar GM [km³/s²]
# prt3("GM of Sun = ", mu_s, " km³ s⁻²")
# ================================================================================================

# ==================================================================================================
# ☀️ CENTRAL STAR PARAMETERS AND METRIC INITIALIZATION
# --------------------------------------------------------------------------------------------------
# This section sets up all physical parameters related to the central massive body (the Sun),
# derives its Schwarzschild and Kerr parameters, and computes the gravitational and
# rotational constants entering the metric tensor used in the integration.
#
# The final part also initializes the state vector for the spacecraft (y0_pt),
# which will later be integrated along the chosen metric.
# ==================================================================================================


# ================================================================================================
# 🌞 GRAVITATIONAL CONSTANTS FOR THE CENTRAL STAR
# ------------------------------------------------------------------------------------------------
# mu — gravitational parameter (G × M)
# Mstar — stellar mass
# rs — Schwarzschild radius (2GM / c²)
# ================================================================================================
const mu = mu_s  # GM of the Sun (in km³/s²)
# To switch to flat Minkowski spacetime for test purposes:
# mu = BF0
prt3("GM of star = ", Float64(mu), " km³ s⁻²")

# Derive stellar mass (kg)
const Mstar = mu / G
prt3("mass of star = ", Float64(Mstar), " kg")

# Schwarzschild radius (in km)
const rs = eps_c2 * BF2 * mu
prt3("Schwarzschild radius of star = ", Float64(rs), " km")
# ================================================================================================


# ================================================================================================
# 🌀 SELECTING THE ACTIVE METRIC
# ------------------------------------------------------------------------------------------------
# The simulation can use either the Kerr (rotating) or Schwarzschild (static) metric.
# The name is retrieved through the Give_Function_Name() mechanism for runtime identification.
# ================================================================================================
const Given_Metric = Give_Function_Name(Kerr_Metric_Polar)
prt3("Given_Metric = ", Given_Metric)
# ================================================================================================


# ================================================================================================
# ⚙️ KERR METRIC INITIALIZATION (ROTATING CENTRAL BODY)
# ------------------------------------------------------------------------------------------------
# When the chosen metric is Kerr, additional rotational parameters are introduced:
#   - Jsun : solar angular momentum (kg·km²/s)
#   - Jstar : same value, reused for consistency
#   - as : scaled angular momentum per unit mass (in km), i.e. a = J / (M c)
#   - event_horizon : radius of the outer Kerr event horizon
# ================================================================================================
if Given_Metric == Kerr_Metric_Polar_name
    prt3("initialization of Kerr metric")

    # Solar angular momentum (from helioseismology)
    # Mauro et al. (2000), ADS 2000ASPC..198..353D
    const Jsun = parse(BigFloat, "2.02e35")  # [kg·km²/s]
    prt3("angular momentum of the Sun = ", Float64(Jsun), " km² kg s⁻¹")

    # For generality, rename to Jstar
    const Jstar = Jsun
    prt3("angular momentum of star = ", Float64(Jstar), " km² kg s⁻¹")

    # Scaled angular momentum (a = J / (M c)), in km
    const as = Jstar / Mstar / c_light
    # To revert to a non-rotating Schwarzschild metric, use:  const as = BF0
    prt3("scaled angular momentum of star = ", Float64(as), " km")

    # Kerr event horizon (outer radius)
    const event_horizon = rs / BF2 + sqrt(rs^2 / BF4 - as^2)
    prt3("event horizon of star = ", Float64(event_horizon), " km")

else
    # For non-rotating metric
    as = BF0
end
# ================================================================================================


# ================================================================================================
# ☀️ SOLAR RADIATION PRESSURE
# ------------------------------------------------------------------------------------------------
# Defines the radiation pressure exerted by sunlight at 1 AU.
# Reference: IAU nominal solar constant = 1361 W/m²
# The resulting pressure (in N/m²) is obtained from:
#     P = (SolarConstant) / (c × 1000)
# since c is expressed in km/s and we convert to m/s.
# ================================================================================================
const solar_constant = parse(BigFloat, "1361.0")  # W/m²
solar_radiation_pressure = solar_constant / c_light / BF1000  # N/m²
prt3("solar_radiation_pressure = ", Float64(solar_radiation_pressure * BF1000 * BF1000), " micro-N/m²")
# ================================================================================================


# ================================================================================================
# 🚀 INITIAL STATE VECTOR
# ------------------------------------------------------------------------------------------------
# The complete system state vector (y0_pt) contains:
#   - 8 variables for the equations of motion
#   - 2 for proper time evolution
#   - 16 for the tetrad components
#
# The following initialization is for the Parker Solar Probe case,
# but can be replaced freely for other spacecraft or test scenarios.
# ================================================================================================
y0_pt = fill(BF0, nsyst)  # initial state vector (contravariant, proper-time parameterized)
# ================================================================================================

# ==================================================================================================
# 🚀 PARKER SOLAR PROBE PHYSICAL PARAMETERS AND INITIAL STATE STRUCTURE
# --------------------------------------------------------------------------------------------------
# This section defines the spacecraft’s physical characteristics (mass, geometry, reflectivity)
# and describes in detail the structure and meaning of the initial state vector (y0_pt).
#
# All dynamic quantities are defined with respect to *proper time τ*,
# and expressed in the Kerr polar metric reference frame.
# ==================================================================================================


# ================================================================================================
# 🛰️ PARKER SOLAR PROBE PHYSICAL CONSTANTS
# ------------------------------------------------------------------------------------------------
# These constants are taken from mission data and are used for
# solar radiation pressure modeling and trajectory dynamics.
#
# References:
#   - NASA Parker Solar Probe documentation
#   - IAU standard solar flux (≈1361–1367 W/m²)
# ================================================================================================

# --- Spacecraft mass (kg)
const mass_PSP = parse(BigFloat, "655.0")

# --- Solar shield (thermal protection system)
const shield_area = parse(BigFloat, "4.0")  # m²
const shield_reflection_coef = parse(BigFloat, "1.8")

# --- Notes:
#     Solar shield always faces the Sun at periapsis,
#     reflection coefficient > 1 accounts for partial diffuse reflection.
#
#     Typical additional spacecraft parameters (optional):
#       Solar panel reflection coefficient = 1.38
#       Solar panel area = 1.6 m²
#
#     Solar flux at 1 AU = 1367 W/m² (slightly variable over the solar cycle).
# ================================================================================================


# ================================================================================================
# 🧩 INITIAL STATE VECTOR (y0_pt)
# ------------------------------------------------------------------------------------------------
# y0_pt is the **initial state vector** expressed with respect to the *proper time τ*
# of the spacecraft.
#
# It contains 26 components organized as follows:
#   [1–8]   → equations of motion (time and spatial coordinates)
#   [9–10]  → proper and integration time tracking
#   [11–26] → tetrad components (4×4 matrix stored sequentially)
#
# Each element is defined in the Kerr polar metric framework.
# ------------------------------------------------------------------------------------------------
# Component definitions:
#   y0_pt[1]  = t₀         # coordinate time (must always be coordinate time)
#   y0_pt[2]  = (dt/dτ)₀   # initial rate of coordinate time wrt proper time
#   y0_pt[3]  = r₀         # initial radial distance
#   y0_pt[4]  = (dr/dτ)₀   # initial radial velocity (proper time)
#   y0_pt[5]  = φ₀         # initial longitude
#   y0_pt[6]  = (dφ/dτ)₀   # initial angular velocity (proper time)
#   y0_pt[7]  = θ₀         # initial colatitude (π/2 for equatorial plane)
#   y0_pt[8]  = (dθ/dτ)₀   # initial polar velocity
#   y0_pt[9]  = τ₀         # proper time (usually 0)
#   y0_pt[10] = evpar₀     # integration time parameter (usually 0)
#
# Remaining 16 elements [11–26] contain the Fermi–Walker tetrad (4×4 identity matrix at t₀).
# ================================================================================================

# Initialize empty 26-element vector
y0_pt = fill(BF0, nsyst)

# ================================================================================================
# 🔗 (Optional) DE430 EPHEMERIDES INTERFACE
# ------------------------------------------------------------------------------------------------
# The DE430 planetary ephemerides can be loaded for real mission conditions using:
#     BODY_GM = init_de430()
# This replaces the simplified Keplerian initialization.
# ================================================================================================
# BODY_GM = init_de430()  # Uncomment to activate DE430 interface
# ================================================================================================
# ==================================================================================================
# 🚀 INITIAL CONDITIONS — PARKER SOLAR PROBE (UT 2025-06-09T00:00:00)
# --------------------------------------------------------------------------------------------------
# This section defines the initial state of the Parker Solar Probe (PSP) in heliocentric coordinates.
#
# The vector `y0_pt` represents the spacecraft’s 4-position and 4-velocity in the Kerr polar metric,
# expressed with respect to its **proper time τ**.
#
# These values correspond to a pure heliocentric orbit (no planetary perturbations)
# and can be replaced with DE430 ephemeris data for real mission initialization.
# ==================================================================================================


# ================================================================================================
# 🕒 INITIAL TIME AND 4-VELOCITY NORMALIZATION
# ------------------------------------------------------------------------------------------------
# y0_pt[1] — coordinate time t₀ (seconds since J2000)
# y0_pt[2] — dt/dτ, the initial derivative of coordinate time wrt proper time
# ================================================================================================
y0_pt[1] = parse(BigFloat,
    "8.026992691847092000000000000000000000000000000000000000000000000000000013e+08"
)  # J2000 epoch (seconds)
y0_pt[2] = parse(BigFloat,
    "1.000000036065583982352343088267053815632478870312984159711699273598484404"
)  # normalized coordinate-time rate
# ================================================================================================


# ================================================================================================
# 🌍 SPATIAL COORDINATES (POLAR FRAME)
# ------------------------------------------------------------------------------------------------
# y0_pt[3] — radial distance (r₀)
# y0_pt[4] — radial velocity (dr/dτ)
# y0_pt[5] — longitude (φ₀)
# y0_pt[6] — angular velocity (dφ/dτ)
# y0_pt[7] — colatitude (θ₀) = π/2 for equatorial orbit
# y0_pt[8] — polar velocity (dθ/dτ)
# ================================================================================================

# --- Radial position: ~60 million km (≈0.4 AU)
y0_pt[3] = parse(BigFloat,
    "6.05555029296586368735956299733688278740147641652955018507419361797892692e+07"
)

# --- Radial velocity (negative: approaching perihelion)
y0_pt[4] = parse(BigFloat,
    "-39.20853470649714613562963656044536896149439334319005675429975900686120154"
)

# --- Longitude (radians)
y0_pt[5] = parse(BigFloat,
    "-0.1051484868287112526396508187748939473128207329184157695124121577662568626"
)

# --- Angular velocity (radians per proper second)
y0_pt[6] = parse(BigFloat,
    "3.269982323577800306566017661024433703896072724496301294340993590072270686e-07"
)

# --- Colatitude (π/2 → equatorial plane)
y0_pt[7] = parse(BigFloat,
    "1.570796111324400805115005611947296536407776016865991878896744220327295114"
)

# --- Polar angular velocity
y0_pt[8] = parse(BigFloat,
    "2.154706430982307446885933453094922768491406204024097910843647327955235975e-07"
)
# ================================================================================================
# ==================================================================================================
# 🧮 CONSISTENCY CHECKS AND FINALIZATION OF INITIAL STATE VECTOR
# --------------------------------------------------------------------------------------------------
# This section converts the initial proper-time derivatives to coordinate-time derivatives,
# verifies the normalization condition of the 4-velocity, and initializes the proper and
# integration time variables.
# ==================================================================================================


# ================================================================================================
# ⚙️ CONVERSION OF VELOCITIES (PROPER TIME → COORDINATE TIME)
# ------------------------------------------------------------------------------------------------
# The Newtonian version of the state vector (y0_pt_Newton) expresses velocities
# with respect to **coordinate time (t)** instead of **proper time (τ)**.
#
# Conversion rule:
#     v_coord_time = v_proper_time / (dt/dτ)
# ================================================================================================
y0_pt_Newton = fill(BF0, nsyst)  # initialize empty vector

for i = 1:4
    y0_pt_Newton[2*i] = y0_pt[2*i] / y0_pt[2]
end
# ================================================================================================


# ================================================================================================
# 🧾 NORMALIZATION CHECK OF INITIAL 4-VELOCITY
# ------------------------------------------------------------------------------------------------
# The GRAPE function `der1_tau_internal` computes dτ/dτ (i.e., normalization of the 4-velocity)
# within the metric framework provided (here Kerr_Metric_Polar).
#
# This quantity must be numerically equal to **1** for the initial vector to be valid.
# If it differs significantly, the initial data may come from a different metric
# than the one used in this simulation.
# ================================================================================================
dtaudevpar_internal = der1_tau_internal(y0_pt, Kerr_Metric_Polar)
prt3("the value at the following line must be equal 1 wrt numerical precision for the initial velocities of the s/c to be physically acceptable:")
prt3("partial derivative of proper time wrt proper time at start (initial vector) = ", dtaudevpar_internal)
# ================================================================================================


# ================================================================================================
# ⏱️ LINK BETWEEN PROPER TIME AND INTEGRATION TIME
# ------------------------------------------------------------------------------------------------
# The external user-defined function `evpar2tau_external()` defines the relationship between:
#   - the *integration time parameter* (evpar)
#   - and the *proper time* (τ)
#
# Depending on the chosen integration mode (`time_option`), this can be:
#   1 → proper time, 2 → function of proper time,
#   3 → coordinate time, 4 → function of coordinate time.
# ================================================================================================
const tau0 = BF0  # initial proper time (starts at zero)
dtaudevpar_external = evpar2tau_external(time_option, evpar0, y0_pt, Given_Metric)
prt3("partial derivative of proper time wrt integration time at start (imposed) = ", dtaudevpar_external)
# ================================================================================================


# ================================================================================================
# 🧩 FINAL ASSIGNMENTS TO THE INITIAL STATE VECTOR
# ------------------------------------------------------------------------------------------------
# These two entries mark the start of the integration:
#   - y0_pt[9]  → proper time τ₀
#   - y0_pt[10] → integration time parameter (evpar₀)
# ================================================================================================
y0_pt[9]  = tau0    # proper time
y0_pt[10] = evpar0  # integration time
# ================================================================================================
# ==================================================================================================
# 🧭 INITIALIZATION OF THE COMOVING TETRAD (FERMI–WALKER TRANSPORT)
# --------------------------------------------------------------------------------------------------
# The tetrad defines a local orthonormal reference frame comoving with the spacecraft.
# Its time-like axis is aligned with the 4-velocity of the probe,
# ensuring a Fermi–Walker transport of orientation throughout the integration.
#
# This process is crucial to correctly simulate on-board measurements
# (e.g., Doppler shifts, accelerometer readings) in the local inertial frame.
# ==================================================================================================


# ================================================================================================
# 🧮 DEFINITION OF LINEARLY INDEPENDENT 4-VECTORS
# ------------------------------------------------------------------------------------------------
# We define 4 linearly independent vectors to serve as the input basis for
# the Gram–Schmidt orthogonalization process.
#
# The first vector (column 1) **must** be proportional to the 4-velocity
# to ensure proper Fermi–Walker transport.
# ================================================================================================

SetLinIndVect = fill(BF0, (4, 4))

# --- First column: proportional to the 4-velocity (u^μ)
SetLinIndVect[1, 1] = y0_pt[2]
SetLinIndVect[2, 1] = y0_pt[4]
SetLinIndVect[3, 1] = y0_pt[6]
SetLinIndVect[4, 1] = y0_pt[8]

# --- Remaining basis vectors: canonical spatial unit directions
SetLinIndVect[1, 2] = BF0 ; SetLinIndVect[2, 2] = BF1 ; SetLinIndVect[3, 2] = BF0 ; SetLinIndVect[4, 2] = BF0
SetLinIndVect[1, 3] = BF0 ; SetLinIndVect[2, 3] = BF0 ; SetLinIndVect[3, 3] = BF1 ; SetLinIndVect[4, 3] = BF0
SetLinIndVect[1, 4] = BF0 ; SetLinIndVect[2, 4] = BF0 ; SetLinIndVect[3, 4] = BF0 ; SetLinIndVect[4, 4] = BF1
# ================================================================================================


# ================================================================================================
# 🧩 COMPUTATION OF THE INITIAL TETRAD (AND ITS INVERSE)
# ------------------------------------------------------------------------------------------------
# The `Vierbein` function orthonormalizes the input basis vectors according to the metric tensor,
# yielding the Fermi–Walker tetrad and its inverse at the initial position.
# ================================================================================================
(tetrad, tetradinv) = Vierbein(SetLinIndVect, Given_Metric, y0_pt)
# ================================================================================================


# ================================================================================================
# 🔗 STORE THE TETRAD INSIDE THE INITIAL STATE VECTOR
# ------------------------------------------------------------------------------------------------
# The tetrad is embedded within the global state vector y0_pt so that
# the integrator can propagate both the motion and the local frame together.
# ================================================================================================
y0_pt = y_from_FWtetrad(y0_pt, tetrad)
# ================================================================================================


# ================================================================================================
# 🧾 VALIDATION OF THE TETRAD ORTHONORMALITY
# ------------------------------------------------------------------------------------------------
# The tetrad is verified using the function `verif_tetrad`, which returns
# the differences (dift, difti) between the computed and expected metric identities.
#
# Both dift and difti should be small (e.g., < 1e-10) to confirm correct orthonormality.
# ================================================================================================
(tetrad, ) = FWtetrad_from_y(y0_pt)  # optional verification
(dift, difti) = verif_tetrad(tetrad, y0_pt, Given_Metric)
prt3("validation of Fermi-Walker (FW) tetrad at start = ",
     Float64(log10(dift)), Float64(log10(difti)), "(log10 of precision)")
# ================================================================================================


# ================================================================================================
# 🚀 COMPUTATION OF INITIAL 4-VELOCITY AND CHOLESKY REPRESENTATION
# ------------------------------------------------------------------------------------------------
# The 4-velocity (in chart coordinates) is extracted from the initial state vector,
# and the metric tensor is evaluated at that position to verify the invariant normalization.
#
# The state is then transformed into its Cholesky-tetrad representation (`ychol0_pt`),
# which is used internally by the GRAPE integrator for numerical stability.
# ================================================================================================
chart_4_velocity = velocity_from_y(y0_pt)
(mt_s, mt_si) = Given_Metric(y0_pt)
invariant0 = mt_s * chart_4_velocity
ychol0_pt = ychol_from_ychart(y0_pt, Given_Metric)
# ================================================================================================


# ================================================================================================
# 📡 DOPPLER REFERENCE FREQUENCY
# ------------------------------------------------------------------------------------------------
# This frequency is used later for Doppler shift simulations.
# Default: 8 GHz (X-band, typical for deep-space communications).
# ================================================================================================
transmit_frequency_Doppler = parse(BigFloat, "8.0e9")
prt3("transmit frequency for Doppler = ", Float64(transmit_frequency_Doppler), "Hertz")
# ================================================================================================

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# ==================================================================================================
# 🧩 MAIN INTEGRATION PHASE
# --------------------------------------------------------------------------------------------------
# This section performs the actual numerical integration of the relativistic trajectory.
# It initializes the symplectic (Cholesky) tetrad, runs the core GRAPE integrator,
# and checks the validity of tetrads and conserved quantities before and after the run.
# ==================================================================================================


# ================================================================================================
# ⏱️ INTEGRATION START — TIME LOGGING
# ------------------------------------------------------------------------------------------------
# Record the computer clock time for performance and traceability.
# ================================================================================================
t_start = now()
prt3("integration start at = ", t_start)
# ================================================================================================

# ================================================================================================
# 🔁 INTEGRATION LOOP WITH DATA COLLECTION  (fixed global scope)
# ================================================================================================
r_plot = Float64[]
phi_plot = Float64[]
tau_plot = Float64[]
time_coord_plot = Float64[]

# define globals you'll mutate in the loop
y = copy(y0_pt)
evpar = evpar0

for k in 1:N_integration_steps
    global y
    global evpar

    y = rki(evpar, evpar + evpar_step,
            y, Given_Metric,
            Eqs_Motion_chol_simple,
            Non_Grav_Grad,
            Non_Grav_Bus)

    evpar += evpar_step

    # collect samples for plots/diagnostics
    push!(r_plot, Float64(y[3]))
    push!(phi_plot, Float64(mod2pi(y[5])))
    push!(tau_plot, Float64(y[9]) / 86400.0)
    push!(time_coord_plot, Float64(y[1]) / 86400.0)
end

# expose final values for the rest of the script
y_end = y
evpar_end = evpar

println("✅ Données d’intégration stockées : $(length(r_plot)) points.")
println("✅ Intégration terminée.")
println("Durée intégrée = ", Float64(evpar_end - evpar0), " s")
println("État final :")
display(y_end)


# ================================================================================================
# 📊 FINALIZATION — OUTPUT SUMMARY (no tetrad validation)
# ================================================================================================
prt3("coordinate time at start and end = ",
     Float64(y0_pt[1]), " s, ",
     Float64(y_end[1]), " s, Δ = ",
     Float64(y_end[1] - y0_pt[1]))

prt3("proper time at start and end = ",
     Float64(y0_pt[9]), " s, ",
     Float64(y_end[9]), " s, Δ = ",
     Float64(y_end[9] - y0_pt[9]))

prt3("integration time at start and end = ",
     Float64(y0_pt[10]), " s, ",
     Float64(y_end[10]), " s, Δ = ",
     Float64(y_end[10] - y0_pt[10]))

dtaudevpar_internal_end = der1_tau_internal(y_end, Given_Metric)
prt3("partial derivative of proper time wrt proper time at end (propagated) = ",
     dtaudevpar_internal_end)

dtaudevpar_external_end = evpar2tau_external(time_option, evpar_end, y_end, Given_Metric)
prt3("partial derivative of proper time wrt integration time at end (imposed) = ",
     dtaudevpar_external_end)
# ================================================================================================

# ================================================================================================
# 🧭 INVARIANT CHECKS — CONSERVATION OF MOMENTA
# ------------------------------------------------------------------------------------------------
# Compute the covariant 4-momentum components at the end of integration
# and compare them with initial values. In a pure Kerr spacetime (no
# non-gravitational forces), the 1st and 3rd components are constants of motion.
# ================================================================================================
chart_4_velocity = velocity_from_y(y_end)
(mt_s, mt_si) = Given_Metric(y_end)  # metric tensor (covariant and contravariant)
invariant = mt_s * chart_4_velocity

prt3("covariant momenta at start and end of integration ==> next 5 lines")
prt3("first and third components are constants of motion for a pure Kerr spacetime with no non-gravitational forces")

for i = 1:4
    prt3(i, invariant[i], invariant0[i])
end
# ================================================================================================


# ================================================================================================
# 📐 NORMALIZATION CHECK IN THE FINAL LOCAL TETRAD FRAME
# ------------------------------------------------------------------------------------------------
# Compute the 4-velocity in the final tetrad frame and check its Minkowski norm.
# This value should remain close to 1 (in units of c_light²) for consistency.
# ================================================================================================
local_velocity = tetrad_vector_from_chart_vector(chart_4_velocity, y_end)
norm2 = Minkowski_norm2(local_velocity)
prt3("control 4-velocity norm in Fermi-Walker tetrad, should be close to 1 (wrt c_light_software^2) = ",
     norm2 / c_light^2)
# ================================================================================================


# ================================================================================================
# 🕒 END OF INTEGRATION — LOGGING
# ------------------------------------------------------------------------------------------------
# Record the computer clock at the end of the run and log elapsed time.
# ================================================================================================
t_end = now()
prt3("computation end = ", t_end)
# ================================================================================================

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# ==================================================================================================
# 🧾 POST-INTEGRATION DIAGNOSTICS AND CONSISTENCY CHECKS
# --------------------------------------------------------------------------------------------------
# After the main integration loop, we compare all time parameters (coordinate, proper, integration),
# output the initial/final state vectors, and verify the derivative relationships between
# proper time and integration time. The elapsed wall-clock time is also reported.
# ==================================================================================================


# ================================================================================================
# 🕰️ COMPARISON OF COORDINATE, PROPER, AND INTEGRATION TIMES
# ------------------------------------------------------------------------------------------------
# These quantities are expected to evolve consistently across the integration.
# Their differences give a first indication of the relativistic effects.
# ================================================================================================
prt3("coordinate time at start and end = ",
     Float64(y0_pt[1]), "seconds, ",
     Float64(y_end[1]), "seconds, ",
     "dif = ", Float64(y_end[1] - y0_pt[1]))

prt3("proper time at start and end = ",
     Float64(y0_pt[9]), "seconds, ",
     Float64(y_end[9]), "seconds, ",
     "dif = ", Float64(y_end[9] - y0_pt[9]))

prt3("integration time at start and end = ",
     Float64(y0_pt[10]), "seconds, ",
     Float64(y_end[10]), "seconds, ",
     "dif = ", Float64(y_end[10] - y0_pt[10]))
# ================================================================================================


# ================================================================================================
# 🧩 FULL STATE VECTOR COMPARISON
# ------------------------------------------------------------------------------------------------
# Display each component of the state vector (position, velocity, tetrad, etc.)
# before and after integration. This provides a detailed diagnostic of numerical
# drift or instability, useful for debugging or validation.
# ================================================================================================
for i = 1:nsyst
    prt3("initial/final state vector (chart) = ", i, y0_pt[i], y_end[i])
end
# ================================================================================================


# ================================================================================================
# 🧮 VALIDATION OF PROPER TIME DERIVATIVES
# ------------------------------------------------------------------------------------------------
# The following quantities verify whether the evolution of proper time (τ)
# remains consistent between internal propagation and external user definitions.
#
#   - `dtaudevpar_internal_end` → computed from metric (should ≈ 1)
#   - `dtaudevpar_external_end` → imposed by user-defined relation
# ================================================================================================
dtaudevpar_internal_end = der1_tau_internal(y_end, Given_Metric)  # internal check
dtaudevpar_external_end = evpar2tau_external(time_option, evpar_end, y_end, Given_Metric)  # external model

prt3("partial derivative of proper time wrt proper time at end (propagated) = ",
     dtaudevpar_internal_end)
prt3("partial derivative of proper time wrt integration time at end (imposed) = ",
     dtaudevpar_external_end)
# ================================================================================================


# ================================================================================================
# ⏱️ PERFORMANCE MEASUREMENT
# ------------------------------------------------------------------------------------------------
# Compute and log the wall-clock time required for the full relativistic integration.
# This provides a direct measure of numerical performance.
# ================================================================================================
prt3("elapsed time for integration = ", t_end - t_start)
# ================================================================================================

# ==================================================================================================
# 📊 PLOTS AND POST-PROCESSING SECTION
# --------------------------------------------------------------------------------------------------
# This section reads the spacecraft ephemerides computed by GRAPE, reconstructs relevant
# physical quantities (position, velocity, Doppler, tetrad accuracy, etc.), and prepares data
# vectors for plotting using PythonPlot.
# ==================================================================================================
using PythonPlot

# --- Rayon en fonction du temps propre ---
figure()
plot(tau_plot, r_plot)
xlabel("Temps propre τ (jours)")
ylabel("Rayon r (km)")
title("Parker Solar Probe — GRAPE NoTetrad")
grid(true)
savefig("orbit_radius_vs_tau_NoTetrad_$(timetag).png")

# --- Projection polaire (r, φ) ---
figure()
plot(r_plot .* cos.(phi_plot), r_plot .* sin.(phi_plot))
xlabel("x (km)")
ylabel("y (km)")
title("Orbital Plane — PSP NoTetrad")
axis("equal")
grid(true)
savefig("orbit_xy_NoTetrad_$(timetag).png")



# ================================================================================================
# 🏁 FINAL STATUS AND CLEANUP
# ------------------------------------------------------------------------------------------------
prt3("timetag = ", timetag)
msg = Call_System("cmd /c echo computation finished")
prt3(msg)
close(runlogfile)
# ================================================================================================
