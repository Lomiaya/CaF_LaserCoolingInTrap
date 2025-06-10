### DEFINE OTHER PARAMETERS FOR THE SIMULATION ###

sim_type = Float64

# σx_initial = 0.2e-6
# σy_initial = 0.2e-6
# σz_initial = 2e-6
# Tx_initial = 40e-6
# Ty_initial = 40e-6
# Tz_initial = 40e-6

Temperature_initial = 70e-6

σx_initial = 0.1e-6
σy_initial = 0.1e-6
σz_initial = 0.5e-6
Tx_initial = Temperature_initial
Ty_initial = Temperature_initial
Tz_initial = Temperature_initial

@everywhere begin
    waist = 1e-6
    Pwr = 20e-3
    lambda = 776e-9
    I0_trap = 2Pwr / (π * waist^2)
end

@everywhere function Is(x, y, z, P, w, lambda)
    r² = x^2 + y^2
    z_R = π * w^2 / lambda  # Rayleigh range
    wz = w * sqrt(1 + (z / z_R)^2)  # Beam waist at z
    return (2 * P / (π * wz^2)) * exp(-2 * r² / wz^2)
end


@everywhere function ∇Is(x, y, z, P, w, lambda)
    I = Is(x, y, z, P, w, lambda)
    r² = x^2 + y^2
    z_R = π * w^2 / lambda  # Rayleigh range
    wz = w * sqrt(1 + (z / z_R)^2)  # Beam waist at z
    d_wz_dz = (z / z_R^2) * (w^2 / wz)  # Derivative of wz with respect to z
    dI_dx = -4 * x * I / wz^2
    dI_dy = -4 * y * I  / wz^2
    dI_dz = 4 * r² * I * d_wz_dz / wz^3 - 2 * I * d_wz_dz / wz
    return (dI_dx, dI_dy, dI_dz)
end

using Serialization
H_ODT_matrix = deserialize("H_tweezer_matrix.jl") * 2π / (2 * ε0 * c) # in unit of s^-1

import MutableNamedTuples: MutableNamedTuple
sim_params = MutableNamedTuple(
    # trap_scalar = sqrt(2I0_trap * 0.0 / (ε0 * c)),
    # H_ODT_matrix = MMatrix{size(H_ODT_matrix)...}(sim_type.(H_ODT_matrix)),
    g = [(x, y, z) -> Is(x, y, z, Pwr, waist, lambda)],
    ∇g = [(x, y, z) -> ∇Is(x, y, z, Pwr, waist, lambda)],
    H₀ = [MMatrix{size(H_ODT_matrix)...}(sim_type.(H_ODT_matrix))],
    
    x_dist = Normal(0, σx_initial),
    y_dist = Normal(0, σy_initial),
    z_dist = Normal(0, σz_initial),
    
    vx_dist = Normal(0, sqrt(kB*Tx_initial/2m)),
    vy_dist = Normal(0, sqrt(kB*Ty_initial/2m)),
    vz_dist = Normal(0, sqrt(kB*Tz_initial/2m)),

    f_z = StructArray(zeros(Complex{sim_type}, 16, 16)),

    dt_diffusion = 1e-7 / (1/Γ)
)