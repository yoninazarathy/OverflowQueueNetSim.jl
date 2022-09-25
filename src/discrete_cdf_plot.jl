"""
Plots a discrete CDF with jumps at x and cumulative values y.
"""
function plot_discrete_cdf(#plt::Plots.Plot{Plots.GRBackend}, 
                            x::AbstractVector, 
                            y::AbstractVector;
                            radius = 0.1, do_scatter = false)
    length(x) != length(y) && error("x and y must be of same length")
    !issorted(x) && error("x must be sorted")
    !issorted(y) && error("y must be sorted")
    
    min_x, max_x = x[1], x[end]
    @show min_x, max_x

    y_adjusted = vcat(0,y[1:end-1])
    y_adjusted_big = vcat(0,y)

    if do_scatter
        plt = scatter(plt, x, y_adjusted,
                        c=:red,
                        markershape=:circle, 
                        markersize = 5,
                        markercolor =nothing, 
                        label = false)
        plt = scatter(plt, x, y,
                        c=:red,
                        markershape=:circle, 
                        markersize = 5,
                        label=false)
    end

    plt = plot([x[1],x[2]],[y[1],y[1]], c=:black,lw=5,label = "Molecule PH approximation" )
    for i in 2:length(x)-1
        plt = plot(plt,[x[i],x[i+1]],[y[i],y[i]], c=:black,lw=3,label =  false)
    end

    for i in 1:length(x)
        plt = plot(plt,[x[i],x[i]],[y_adjusted_big[i],y_adjusted_big[i+1]],c=:black,label=false, linestyle = :dash)
    end
    
    return plt
end