using Distributed
# procs_to_use = 5
# if nprocs() <= procs_to_use
#     addprocs(procs_to_use-nprocs())
# end
@everywhere using QuantumStates, OpticalBlochEquations, DifferentialEquations, UnitsToValue, StructArrays, StaticArrays, Parameters
@everywhere import Distributions: Normal, Geometric, Exponential
@everywhere using ProgressMeter, Plots
@everywhere include("helper_functions.jl")
@everywhere include("define_molecular_structure.jl")
@everywhere include("define_sim_params.jl")
@everywhere include("define_prob.jl")
@everywhere include("compute_size_temperature.jl")

using Serialization
using Printf

# xpoints = [0.6, 0.8, 1.0] # in mW
# ypoints = [0.3, 0.4, 0.5] # in mW
xpoints = [0.2, 0.6, 1.0]
# ypoints = [0.05, 0.15, 0.25, 0.35, 0.45]
# ypoints = [0.15, 0.3, 0.45, 0.6, 0.75]
ypoints = [0.1, 0.2, 0.25]

points = [(i, j) for i in xpoints, j in ypoints]

for element in points
    prob.p.sats = [element[1]*total_sat, element[2]*total_sat]
    prob_diffusion.p.sats = [element[1]*total_sat, element[2]*total_sat]
    global s1 = element[1] * total_sat
    global s2 = element[2] * total_sat

    n_trajectories1 = 50
    n_trajectories2 = 50
    n_times = 10

    n_trajectories_diffusion = 50000 # # of particles ran for the diffusion.
    diffusion_t_end = 0e-6
    diffusion_τ_total = 6e-6

    try
        (sols_no_diffusion, sols_with_diffusion, diffusion, diffusion_error, diffusion_over_time) = 
            compute_trajectories_with_diffusion(
            prob, prob_func!, prob, prob_func_diffusion!, n_trajectories1, n_trajectories2, n_trajectories_diffusion, n_times, diffusion_t_end, diffusion_τ_total
        )
        Ts = T_vs_time(sols_with_diffusion)
        display(Ts[end] .* 1e6)
        xs_dfshb, ys_dfshb = survival_rate_curve(sols_with_diffusion)
        display(ys_dfshb[end])
        display(sols_with_diffusion[10].prob.p.n_scatters / (-log(ys_dfshb[end])))

        filename = @sprintf("results/cooling_sim_d%.2fd%.2fd%.2fb%.2fs%.2fs%.2ft%.2f.jl",
                            detuning,
                            δ1,
                            δ2,
                            beam_radius / 1e-3,
                            s1 / total_sat,
                            s2 / total_sat,
                            Temperature_initial / 1e-6)
        serialize(filename, sols_with_diffusion)

        open("results/output.txt", "a") do io
            println(io, "$filename")
            var_name = "Temperature (μK)"
            formatted_value = @sprintf("%.3g", Ts[end] .* 1e6)
            println(io, "$var_name: $formatted_value")
            var_name = "Scattering rate (count, 3ms)"
            formatted_value = @sprintf("%.3g", sols_with_diffusion[10].prob.p.n_scatters)
            println(io, "$var_name: $formatted_value")
            var_name = "p_survival (count)"
            formatted_value = @sprintf("%.3g", ys_dfshb[end])
            println(io, "$var_name: $formatted_value")
            var_name = "Scattering rate over log p_survival (count)"
            formatted_value = @sprintf("%.3g", sols_with_diffusion[10].prob.p.n_scatters / (-log(ys_dfshb[end])))
            println(io, "$var_name: $formatted_value")
            var_name = "Diffusion constant (e-3)"
            formatted_value = @sprintf("%.3g", diffusion / 1e-3)
            println(io, "$var_name: $formatted_value")
        end
    catch
        continue
    end
end