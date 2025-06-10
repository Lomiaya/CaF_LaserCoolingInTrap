using Serialization
states = deserialize("states_cooling_sim_tweezer.jl")
ground_states = states[1:12]
excited_states = states[13:16]

d = zeros(ComplexF64, 16, 16, 3)
d[1:12, 13:16, :] .= tdms_between_states(ground_states, excited_states)
d[13:16, 1:12, :] .= tdms_between_states(excited_states, ground_states)

# Define constants for the laser cooling transition
@everywhere begin
    @consts begin
        λ = 606e-9
        Γ = 2π * 8.3e6
        m = @with_unit 59 "u"
        k = 2π / λ
    end
end


@everywhere import LoopVectorization: @turbo
@everywhere function add_single_term_dψ!(dψ, ψ, g, H, r, t)
    g_factor = g(r[1] / k,r[2] / k,r[3] / k)
    @turbo for i ∈ 1:16
        dψ_i_re = zero(eltype(dψ.re))
        dψ_i_im = zero(eltype(dψ.im))
        for j ∈ 1:16
            ψ_i_re = ψ.re[j]
            ψ_i_im = ψ.im[j]
            
            H_re = g_factor * H[i,j]
            H_im = 0.0
            
            dψ_i_re += ψ_i_re * H_re - ψ_i_im * H_im
            dψ_i_im += ψ_i_re * H_im + ψ_i_im * H_re
            
        end
        dψ.re[i] += dψ_i_im / Γ
        dψ.im[i] -= dψ_i_re / Γ
    end
    return nothing
end

@everywhere function add_terms_dψ!(dψ, ψ, p, r, t)
    for k ∈ 1:length(p.sim_params.H₀)
        add_single_term_dψ!(dψ, ψ, p.sim_params.g[k], p.sim_params.H₀[k], r, t)
    end
    return nothing
end
