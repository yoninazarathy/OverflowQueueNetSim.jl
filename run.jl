using Pkg
Pkg.activate(".")

using Distributions
using DataStructures
using Parameters
using StatsBase
using Plots
using Random
using LinearAlgebra
using Calculus
using Accessors

include("src/util.jl")
include("src/sim_engine.jl")
include("src/network_def.jl")
include("src/network_action.jl")
include("src/random.jl")
    include("src/approx/ph_molecule_approx.jl")
    include("src/approx/sim_mc_molecule_approx.jl")
    include("src/examples/example1.jl")

save_figs = true
Random.seed!(0)

function do_k(all_k; example = example1, 
                     max_time = 10e5)
    example = @set example.K=all_k*ones(Int,length(example.K))
    sojourn_times, _, _ = do_experiment_traj(example; max_time = max_time, warm_up_time = 10^4)
    return ecdf(sojourn_times/all_k), sojourn_times
end

approx_dist_ecdf = mc_sim_approx_dist(example1, λ1)

scaled_sojourn_ecdf_50, sojourn_times_50 = do_k(50)
scaled_sojourn_ecdf_100, sojourn_times_100 = do_k(100)
scaled_sojourn_ecdf_500, sojourn_times_500 = do_k(500)

top_x = 5
x_grid = -0.01:(top_x/2000):top_x

plot(x_grid, scaled_sojourn_ecdf_50.(x_grid), 
            label = "Buffer Scaling K=50",
            ylim=(0,1),lw=2, c=:green)

plot!(x_grid, scaled_sojourn_ecdf_100.(x_grid), 
            label = "Buffer Scaling K=100",
            ylim=(0,1),lw=2, c=:blue)

plot!(x_grid, scaled_sojourn_ecdf_500.(x_grid), 
            label = "Buffer Scaling K=500",
            ylim=(0,1),lw=2, c=:black)

plot!(x_grid, approx_dist_ecdf.(x_grid), 
            label = "Molecule Approximation",
            legend=:bottomright, c=:red,lw=1,
            xlabel="Sojourn time normalized by K",
            ylabel="Cumulative Distribution")

save_figs && savefig("figs/SojournTimeCDFs.pdf")

top_x = quantile(sojourn_times_50,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times_50, bins =3000,normed=true,
        xlim=(0,top_x),lw=2, c=:black,ylim=(0,0.02),
        xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 50")
save_figs && savefig("figs/SojournTimesHist50.pdf")

top_x = quantile(sojourn_times_100,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times_100, bins =3000,normed=true,
        xlim=(0,top_x),lw=2, c=:black,ylim=(0,0.02),
        xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 100")
save_figs && savefig("figs/SojournTimesHist100.pdf")

top_x = quantile(sojourn_times_500,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times_500,
        xlim=(0,top_x), c=:black, bins = 3000,
        normed=true,ylim=(0,0.02),
        lw = 2,xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 500")
save_figs && savefig("figs/SojournTimesHist500.pdf")
