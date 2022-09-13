using Pkg
Pkg.activate(".")

using Distributions
using DataStructures
using Parameters
using StatsBase
using Plots
using Random
using LinearAlgebra

include("util.jl")
include("sim_engine.jl")
include("network_def.jl")
include("network_action.jl")
include("examples.jl")
include("random.jl")
include("limiting_approx.jl")
include("markov_reward.jl")

outer_α = 1.5
p = OverflowNetworkParameters(  K = ones(Int,5), 
                                α=outer_α*[0.5, 0.1, 0.2, 0.2, 0.4],
                                μ=[1,1,1,1,1],
                                P = [0.0 0.1 0.3 0.1 0.0; 
                                    0.4 0.0 0.4 0.1 0.0; 
                                    0.2 0.1 0.0 0.1 0.2;
                                    0.3 0.1 0.1 0.0 0.0;
                                    0.0 0.2 0.3 0.1 0.0],
                                Q = [0.0 0.1 0.2 0.2 0.0; 
                                    0.2 0.0 0.2 0.3 0.0; 
                                    0.0 0.5 0.0 0.1 0.1;
                                    0.0 0.2 0.2 0.0 0.0;
                                    0.2 0.0 0.1 0.1 0.0],                                    
                                α_scv = 1.0*ones(5),
                                μ_scv = 1.0*ones(5))

# make_init_dist(p)

dist_numeric = make_dist(p)
dist_simulation = sim_approx_dist(p)

[(numeric = dist_numeric(k), sim = dist_simulation(k)) for k in 0:10]

# x_grid = 0.01:0.01:5
# plot(x_grid,dist_numeric.(x_grid))