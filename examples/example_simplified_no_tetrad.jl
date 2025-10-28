#!/usr/bin/env julia
# ============================================================
#  GRAPE Example — Parker Solar Probe (Simplified Test)
#  Avec barre de progression et pseudo-Doppler
#  Compatible avec module GRAPE_P (GRAPE_core.jl)
# ============================================================

# --- 1️⃣ Chargement du moteur physique (GRAPE)
include("../src/GRAPE_core.jl")
using .GRAPE_P
using Dates, Printf, LinearAlgebra, PythonPlot, ProgressMeter, DelimitedFiles

println("🚀 Chargement du module GRAPE_P terminé.")

# --- 2️⃣ Préparation du dossier de sortie
mkpath("output")

# --- 3️⃣ Définition des constantes physiques (unités SI)
const c_light = BigFloat(2.99792458e8)           # vitesse de la lumière (m/s)
const mu      = BigFloat(1.32712440018e20)       # GM du Soleil (m^3/s^2)
const rs      = BigFloat(2 * mu / c_light^2)     # rayon de Schwarzschild du Soleil
const as      = BigFloat(0.0)                    # paramètre de rotation (ici nul)
const AU      = BigFloat(1.495978707e11)         # unité astronomique (m)

println("✅ Constantes physiques initialisées.")

# --- 4️⃣ Paramètres d'intégration
const evpar_step = BigFloat(5.0)                 # pas de temps d’intégration (s)
const N_integration_steps = 1200                 # nombre total de pas
const time_option = 1                            # intégration au temps propre
const order_parder_metric = 4                    # ordre dérivée du tenseur métrique
const itmax_rki = 6                              # itérations méthode implicite
const epsmin_CSN = BigFloat(1e-6)
const epsabs_CSN = BigFloat(1e-6)
const eps_tetrad = BigFloat(1e-10)

println("✅ Paramètres d’intégration définis.")
println("   - Pas de temps : $(Float64(evpar_step)) s")
println("   - Nombre d’étapes : $N_integration_steps")

# --- 5️⃣ Conditions initiales : Parker Solar Probe au périhélie
r0 = BigFloat(9.86 * 6.957e8)                    # rayon initial (~9.86 R☉)
θ0 = BigFloat(pi/2)                              # équatorial
ϕ0 = BigFloat(0.0)
v0_tan = sqrt(mu / r0) * BigFloat(0.6)           # 60 % vitesse circulaire

println("🌞 Initialisation de la sonde Parker Solar Probe...")
println("   Rayon initial      : $(Float64(r0/6.957e8)) R☉")
println("   Vitesse tangente   : $(Float64(v0_tan/1000)) km/s")

# Le vecteur d’état y a 26 composantes (voir module GRAPE_P)
y0 = fill(BigFloat(0), 26)
# coordonnées (t, r, φ, θ)
y0[1] = 0.0
y0[3] = r0
y0[5] = ϕ0
y0[7] = θ0
# vitesses correspondantes
y0[2] = 1.0
y0[4] = 0.0
y0[6] = v0_tan
y0[8] = 0.0
# temps propres et d’intégration
y0[9]  = 0.0
y0[10] = 0.0
# tetrad initial (identité)
for i in 1:4
    for j in 1:4
        y0[4*i+j+6] = i == j ? BigFloat(1) : BigFloat(0)
    end
end

println("✅ État initial défini (format GRAPE).")

# --- 6️⃣ Conversion dans le repère de Cholesky
ychol0 = GRAPE_P.ychol_from_ychart(y0, GRAPE_P.Kerr_Metric_Polar)

# --- 7️⃣ Intégration principale avec barre de progression
println("\n🧮 Début de l’intégration relativiste...")
timetag = Dates.format(now(), "yyyymmdd_HHMMSS")

# Progress bar setup
p = Progress(N_integration_steps; dt=0.5,
             desc="Intégration en cours",
             barlen=40, color=:cyan)

# Ouverture du fichier d’éphémérides
open("SCephemeris_$(timetag).txt", "w") do f
    ychol = deepcopy(ychol0)
    evpar = BigFloat(0.0)
    nsyst = length(ychol)

    # --- pour stocker le pseudo-Doppler à chaque pas
    doppler_ratio = zeros(Float64, N_integration_steps)
    time_values = zeros(Float64, N_integration_steps)

    for i in 1:N_integration_steps
        ychol = GRAPE_P.rki(evpar, evpar + evpar_step, ychol,
                            GRAPE_P.Kerr_Metric_Polar,
                            GRAPE_P.Eqs_Motion_chol,
                            GRAPE_P.Non_Grav_Grad,
                            GRAPE_P.Non_Grav_Bus)
        evpar += evpar_step

        # Écriture de l’état
        for j in 1:nsyst
            write(f, "$(j) $(ychol[j])\n")
        end

        # ⚙️ Calcul du pseudo-Doppler (variable décrite dans GRAPE_core)
        # ychol[2] = dt/dτ (coord. time / proper time)
        # ychol[9] = dτ/dσ (proper / integration time)
        # ychol[10] = dσ (integration time)
        dt_dσ = Float64(ychol[2])   # dérivée du temps coordonné
        dτ_dσ = Float64(ychol[9])   # dérivée du temps propre
        doppler_ratio[i] = dt_dσ / dτ_dσ
        time_values[i] = Float64(evpar)

        next!(p)  # mise à jour de la barre
    end

    finish!(p)
    println("✅ Intégration terminée !")

    # --- Sauvegarde du pseudo-Doppler pour analyse
    doppler_data = hcat(time_values, doppler_ratio)
    writedlm("output/pseudo_doppler_$(timetag).txt", doppler_data)
end

# --- 8️⃣ Lecture des résultats
filename = "SCephemeris_$(timetag).txt"
println("📂 Lecture du fichier de sortie : $filename")
data = readdlm(filename)

# Reshape : le fichier a 26 lignes par pas de temps
idx = reshape(data[:, 1], (26, :))
vals = reshape(data[:, 2], (26, :))

# Extraire les coordonnées polaires et calculer les positions cartésiennes
r_vals = vals[3, :]
ϕ_vals = vals[5, :]
x_vals = r_vals .* cos.(ϕ_vals)
y_vals = r_vals .* sin.(ϕ_vals)

# --- 9️⃣ Affichage graphique
println("🪐 Génération des graphiques...")

# 9a — Trajectoire orbitale
figure(figsize=(6,6))
plot(x_vals ./ (6.957e8), y_vals ./ (6.957e8),
     lw=1.8, color="tab:blue", label="Trajectoire PSP")
xlabel("x / R☉")
ylabel("y / R☉")
title("Parker Solar Probe — Orbite relativiste (GRAPE simplifié)")
grid(true)
axis("equal")
legend()
out_orbit = "output/PSP_orbit_simple.png"
savefig(out_orbit, dpi=200)
println("✅ Graphique d’orbite sauvegardé : $out_orbit")

# 9b — Pseudo-Doppler (affecte la fréquence reçue)
doppler_data = readdlm("output/pseudo_doppler_$(timetag).txt")
time_vals = doppler_data[:,1]
doppler_vals = doppler_data[:,2]

figure(figsize=(6,4))
plot(time_vals ./ 3600, doppler_vals,
     lw=1.8, color="tab:red", label="dt/dτ ratio")
xlabel("Temps d’intégration (heures)")
ylabel("dt/dτ (pseudo-Doppler)")
title("Variation temporelle du pseudo-Doppler (relativiste)")
grid(true)
legend()
out_doppler = "output/PSP_pseudo_doppler.png"
savefig(out_doppler, dpi=200)
println("✅ Graphique Doppler sauvegardé : $out_doppler")

# --- 🔟 Résumé final
println("\n🎯 Simulation complétée avec succès !")
println("   - Fichier d’éphémérides : $filename")
println("   - Image d’orbite         : $out_orbit")
println("   - Image Doppler          : $out_doppler")
println("   - Nombre d’étapes        : $N_integration_steps")
println("   - Pas de temps           : $(Float64(evpar_step)) s")
println("\nFin du test Parker Solar Probe ✅")
