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


