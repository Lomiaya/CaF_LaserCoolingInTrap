@everywhere import LoopVectorization: @turbo

@everywhere function operator_matrix_expectation_complex(O, state)
    O_re = 0.0
    O_im = 0.0
    @turbo for i ∈ 1:16
        re_i = state.re[i]
        im_i = state.im[i]
        for j ∈ 1:16
            re_j = state.re[j]
            im_j = state.im[j]
            cicj_re = re_i * re_j + im_i * im_j # real part of ci* * cj
            cicj_im = re_i * im_j - im_i * re_j
            O_re += O[i,j] * cicj_re - 0.0 * cicj_im
            O_im += O[i,j] * cicj_im + 0.0 * cicj_re
        end
    end
    return (O_re, O_im)
end

@everywhere function add_force_from_hamiltonian!(F, ψ, sim_params, r)
    for i ∈ 1:length(sim_params.H₀)
        ∇g = sim_params.∇g[i](r[1] / k, r[2] / k, r[3] / k)
        H = sim_params.H₀[i]
        E = operator_matrix_expectation_complex(H, ψ)[1] / (k * Γ)
        F[1] -= ∇g[1] * E
        F[2] -= ∇g[2] * E
        F[3] -= ∇g[3] * E
    end
    return nothing
end
