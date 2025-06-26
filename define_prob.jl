### SIMULATION PROBLEMS DEFINITION: 1+2 BLUE MOT ###
include("cooling_params.jl")
include("add_force_from_hamiltonian.jl")

# PROBLEM TO CALCULATE TRAJECTORIES #
t_start = 0.0
t_end   = 3e-3
t_span  = (t_start, t_end) ./ (1/Γ)

p = initialize_prob(sim_type, energies, freqs, sats, pols, beam_radius, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p!, add_terms_dψ!, add_force_from_hamiltonian!)

@everywhere function diffusion_kick(integrator)
    p = integrator.p
    n_states = p.n_states
    n_excited = p.n_excited
    @inbounds @fastmath for i ∈ 1:3
        kick = sqrt(2p.diffusion_constant[i] * p.sim_params.dt_diffusion) / p.m
        integrator.u[2n_states + n_excited + 3 + i] += rand((-1,1)) * kick
    end  
    return nothing
end

diffusion_times = (0:1e-6:t_end) ./ (1/Γ)
cb_diffusion = PresetTimeCallback(diffusion_times, diffusion_kick, save_positions=(false,false))
cb1 = DiscreteCallback(condition_discrete, stochastic_collapse_nokick!, save_positions=(false,false))
cb2 = DiscreteCallback(terminate_condition, terminate!)
cbs = CallbackSet(cb1, cb2, cb_diffusion)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=1000000000, callback=cbs)
prob = ODEProblem(ψ_fast!, p.u0, sim_type.(t_span), p; kwargs...) # this have forces and stuff, have diffusion included, but not the one used to calculate the diffusion constant

# PROBLEM TO COMPUTE DIFFUSION #
p_diffusion = initialize_prob(sim_type, energies, freqs, sats, pols, Inf, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p_diffusion!, add_terms_dψ!, add_force_from_hamiltonian!)
cbs_diffusion = CallbackSet(cb1, cb_diffusion)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=1000000000, callback=cbs_diffusion)
prob_diffusion = ODEProblem(ψ_fast_ballistic!, p_diffusion.u0, sim_type.(t_span), p_diffusion; kwargs...)

function prob_func_CaF!(prob)
    prob.p.n_scatters = 0.
    prob.p.last_decay_time = 0.
    prob.p.time_to_decay = rand(Exponential(1))
    prob.p.sim_params.photon_budget = 1000000
    update_initial_position!(prob, sample_position(prob.p.sim_params) ./ (1/k))
    update_initial_velocity!(prob, sample_velocity(prob.p.sim_params) ./ (Γ/k))
    update_phases!(prob)
    return nothing
end

function prob_func_diffusion_CaF!(prob)
    prob_func_CaF!(prob)
    return nothing
end