function plot_discrete_cdf(plt::Plots.Plot{Plots.GRBackend}, x::Vector{Int}, y::Vector{Float64})
    plt = plot(plt,rand(20))
    return plt
end