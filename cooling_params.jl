# # DEFINE STATES #
# energy_offset = (2π / Γ) * energy(states[13])
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +18.0
# detuning = +40.0
# detuning = +32.0
δ1 = +0.00 # F = 1-
δ2 = +0.00 # F = 2

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[1]) + Δ1 # F = 1-
f2 = energy(states[end]) - energy(states[10]) + Δ2 # F = 2

freqs = [f1, f2] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
# beam_radius = 1e-3
beam_radius = 2e-3
Isat = π*h*c*Γ/(3λ^3)
# P = 1e-3
P = 1e-3
I = 2P / (π * beam_radius^2)

total_sat = I / Isat
s1 = 0.7total_sat # F = 1-, in mW
s2 = 0.35total_sat # F = 2, in mW
# s1 + s2 = I

sats = [s1, s2]

(s1, s2) |> display

# DEFINE POLARIZATIONS #
pols = [σ⁻, σ⁺]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)
    return nothing
end

@everywhere function update_p_diffusion!(p, r, t)
    return nothing
end
