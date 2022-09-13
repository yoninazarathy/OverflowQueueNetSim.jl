using Pkg
Pkg.activate(".")

using Distributions
using DataStructures
using Parameters
using StatsBase
using Plots
# using StatsPlots
using Random
using LinearAlgebra
using Calculus

include("src/util.jl")
include("src/sim_engine.jl")
include("src/network_def.jl")
include("src/network_action.jl")
include("src/examples.jl")
include("src/random.jl")
# include("limiting_approx.jl")
include("src/markov_reward.jl")

# anim = Animation()
# for all_k in [10,10,10,50,100,200,500,1000,2000,5000,5000]

Random.seed!(0)

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
α_scv = 3.0*ones(5),
μ_scv = 0.7*ones(5))


function do_k(all_k; scaled = true, max_time = 10e5)

    outer_α = 1.5
    p = OverflowNetworkParameters(  K = all_k*ones(Int,5), 
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




    sojourn_times, queues, rates = do_experiment_traj(p;max_time=max_time,warm_up_time = 10^4)
    # @show rates
    # zero_sojourn_prob = round(mean(iszero.(sojourn_times)),digits=2)
    # @show zero_sojourn_prob
    # queue_occupancy =  round.(queues/all_k,digits=2)
    # @show queue_occupancy
    # println(summarystats(sojourn_times))
    return scaled ? ecdf(sojourn_times/all_k) : sojourn_times
end

approx_dist_ecdf = sim_approx_dist(p)


top_x = 5#quantile(sojourn_times,0.99)
x_grid = -0.01:(top_x/2000):top_x
    
sojourn_ecdf_30 = do_k(30)
plot(x_grid, sojourn_ecdf_30.(x_grid), label = "Buffer Scaling K=30",
            ylim=(0,1),lw=2, c=:green)

sojourn_ecdf_300 = do_k(300)
plot!(x_grid, sojourn_ecdf_300.(x_grid), label = "Buffer Scaling K=300",
            ylim=(0,1),lw=2, c=:blue)

sojourn_ecdf_3000 = do_k(3000)
plot!(x_grid, sojourn_ecdf_3000.(x_grid), label = "Buffer Scaling K=3000",
            ylim=(0,1),lw=2, c=:black)

plot!(x_grid, approx_dist_ecdf.(x_grid), 
            label = "Molecule Approximation",
            legend=:bottomright, c=:red,lw=1,
            xlabel="Sojourn time normalized by K",
            ylabel="Cumulative Distribution")

savefig("figs/SojournTimeApprox.pdf")

# frame(anim)
# gif(anim, "varying_k.gif", fps = 2)

sojourn_times = do_k(50; scaled = false,max_time =10e5)
top_x = quantile(sojourn_times,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times, bins =3000,normed=true,
        xlim=(0,top_x),lw=2, c=:black,ylim=(0,0.02),
        xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 50")
savefig("figs/SojournTimesHist50.pdf")

sojourn_times = do_k(100; scaled = false,max_time =10e5)
top_x = quantile(sojourn_times,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times, bins =3000,normed=true,
        xlim=(0,top_x),lw=2, c=:black,ylim=(0,0.02),
        xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 100")
savefig("figs/SojournTimesHist100.pdf")


sojourn_times = do_k(500; scaled = false)
top_x = quantile(sojourn_times,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times,
        xlim=(0,top_x), c=:black, bins = 3000,
        normed=true,ylim=(0,0.02),
        lw = 2,xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 500")
savefig("figs/SojournTimesHist500.pdf")

sojourn_times = do_k(1000; scaled = false)
top_x = quantile(sojourn_times,0.99)
x_grid = 0.1:(top_x/100):top_x
histogram(sojourn_times,
        xlim=(0,top_x), c=:black, bins = 3000,
        normed=true,ylim=(0,0.02),
        lw = 2,xlabel="Sojourn time", ylabel="Frequency",label=false,title="Buffer sizes = 1000")
savefig("figs/SojournTimesHist1000.pdf")