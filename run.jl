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
using PyCall

include("src/setup_python.jl")
setup_python() #create solve_overflow_traffic_equation()

include("src/util.jl")
include("src/sim_engine.jl")
include("src/network_def.jl")
include("src/network_action.jl")
include("src/random.jl")
include("src/discrete_cdf_plot.jl")
    include("src/approx/ph_molecule_approx.jl")
    include("src/approx/sim_mc_molecule_approx.jl")
    include("src/examples/example.jl")


function do_sim_k(all_k, example; max_time = 2e6, seed = 0)
        Random.seed!(seed)
        example = @set example.K = all_k*ones(Int,length(example.K))
        sojourn_times, _, _ = do_experiment_traj(example; max_time = max_time, warm_up_time = 10^5)
        return ecdf(sojourn_times/all_k), sojourn_times/all_k
end

function create_cdfs(   example::OverflowNetworkParameters, 
                        λsol::Vector{Float64}, 
                        title::String; 
                        save_figs = true,
                        legend = :bottomright)

        approx_dist_ecdf_ph_type = make_dist(example, λsol)
        # approx_dist_ecdf_mc_sim = mc_sim_approx_dist(example, λsol)

        scaled_sojourn_ecdf_50, sojourn_times_50 = do_sim_k(50, example)
        scaled_sojourn_ecdf_100, sojourn_times_100 = do_sim_k(100, example)
        scaled_sojourn_ecdf_500, sojourn_times_500 = do_sim_k(500, example)
        # scaled_sojourn_ecdf_2000, sojourn_times_2000 = do_sim_k(2000, example)

        top_x = 5
        x_grid = -0.01:(top_x/2000):top_x

        plt = plot_discrete_cdf(0:top_x, [approx_dist_ecdf_ph_type(k) for k in 0:top_x])

        plt = plot(plt,x_grid, scaled_sojourn_ecdf_50.(x_grid), 
                label = "Buffer Scaling K=50",
                ylim=(0,1),lw=2, c=:purple,
                xlabel="Sojourn time normalized by K",
                ylabel="Cumulative Distribution",
                title = title)

        plt = plot(plt, x_grid, scaled_sojourn_ecdf_100.(x_grid), 
                label = "Buffer Scaling K=100",
                ylim=(0,1),lw=2, c=:blue)

        plt = plot(plt, x_grid, scaled_sojourn_ecdf_500.(x_grid), 
                label = "Buffer Scaling K=500",
                ylim=(0,1),lw=2,
                legend=legend, c=:green)

        # plt = plot(plt, x_grid, scaled_sojourn_ecdf_2000.(x_grid), 
        #         label = "Buffer Scaling K=2000",
        #         ylim=(0,1),lw=1.5,
        #         legend=:bottomright, c=:red)


        # plt = plot(plt, x_grid, approx_dist_ecdf_ph_type.(x_grid), 
        #             label = "Molecule Approximation PH",
        #             legend=:bottomright, c=:red,lw=1,
        #             xlabel="Sojourn time normalized by K",
        #             ylabel="Cumulative Distribution")

        # plt = plot(plt, x_grid, approx_dist_ecdf_mc_sim.(x_grid), 
        #                 label = "Molecule Approximation MC Sim",
        #                 legend=:bottomright, c=:purple,lw=1)
                        
        save_figs && savefig(plt,"figs/$(example.name)_SojournTimeCDFs.pdf")
        return nothing
end

function create_sojourn_densities(example::OverflowNetworkParameters; save_figs = true)
        _, sojourn_times_50 = do_sim_k(50, example)
        _, sojourn_times_100 = do_sim_k(100, example)
        _, sojourn_times_500 = do_sim_k(500, example)
        _, sojourn_times_2000 = do_sim_k(2000, example)

        top_x = 5#quantile(sojourn_times_50,0.99)
        x_grid = 0.1:(top_x/100):top_x
        histogram(sojourn_times_50, bins = 3000, normed=true,
                xlim=(0,top_x),lw=2, c=:black, ylim=(0,10),
                xlabel="Sojourn time normalized by K", ylabel="Density",label=false,title="Buffer sizes = 50")
        save_figs && savefig("figs/$(example.name)_SojournTimesHist50.pdf")

        top_x = 5#quantile(sojourn_times_100,0.99)
        x_grid = 0.1:(top_x/100):top_x
        histogram(sojourn_times_100, bins =3000,normed=true,
                xlim=(0,top_x),lw=2, c=:black, ylim=(0, 10),
                xlabel="Sojourn time normalized by K", ylabel="Density",label=false,title="Buffer sizes = 100")
        save_figs && savefig("figs/$(example.name)_SojournTimesHist100.pdf")

        top_x = 5#quantile(sojourn_times_500,0.99)
        x_grid = 0.1:(top_x/100):top_x
        histogram(sojourn_times_500,
                xlim=(0,top_x), c=:black, bins = 3000,
                normed=true,ylim=(0,10),
                lw = 2,xlabel="Sojourn time normalized by K", ylabel="Density",label=false,title="Buffer sizes = 500")
        save_figs && savefig("figs/$(example.name)_SojournTimesHist500.pdf")

        top_x = 5#quantile(sojourn_times_2000,0.99)
        x_grid = 0.1:(top_x/100):top_x
        histogram(sojourn_times_2000,
                xlim=(0,top_x), c=:black, bins = 3000,
                normed=true,ylim=(0,10),
                lw = 2,xlabel="Sojourn time normalized by K", ylabel="Density",label=false,title="Buffer sizes = 2000")
        save_figs && savefig("figs/$(example.name)_SojournTimesHist2000.pdf")

        return nothing
end

println("Running simulations for CDFs of example 1")
create_cdfs(example1, λ1, "$(length(example1)) buffer network with $nfull1 overflow")

println("Running simulations for CDFs of example 2")
create_cdfs(example2, λ2, "$(length(example2)) buffer network with $nfull2 overflows",legend=false)

println("Running simulations for CDFs of example 3")
create_cdfs(example3, λ3, "$(length(example3)) buffer network with $nfull3 overflows",legend=false)

println("Running simulations for CDFs of example 4")
create_cdfs(example4, λ4, "$(length(example4)) buffer network with $nfull4 overflows",legend=false)


println("Running simulations for sojourn densities of example 2")
create_sojourn_densities(example2)