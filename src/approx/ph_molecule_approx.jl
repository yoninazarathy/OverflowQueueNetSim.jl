function make_buffer_order(p::OverflowNetworkParameters, λ::Vector{Float64})
    n = p.n
    full_buffers = (1:n)[λ .> p.μ]
    s = length(full_buffers)
    s_bar = n - s
    buffer_order = vcat(full_buffers, setdiff(1:n,full_buffers))
    return n, s, s_bar, full_buffers, buffer_order
end

function make_absorbtion_matrix(p::OverflowNetworkParameters, λ::Vector{Float64})
    n, s, s_bar, full_buffers, buffer_order = make_buffer_order(p, λ)

    α = p.α[buffer_order]
    P = p.P[buffer_order, buffer_order]
    Q = p.Q[buffer_order, buffer_order]
    γ = [p.μ[i]/λ[i] for i in buffer_order]
    pBar, qBar = 1 .- P*ones(n), 1 .- Q*ones(n)

    C = [Q[1:s,:] .* (1 .- γ[1:s]);
         P[s+1:end, :]]
    
    B = [((1 .- γ[1:s]) .* qBar[1:s]) Diagonal(γ[1:s]); 
         pBar[s+1:end] zeros(s_bar,s)]

    A = inv(Matrix(I-C))*B
    return A
end

function make_dist(p::OverflowNetworkParameters, λ::Vector{Float64})
    n, s, s_bar, full_buffers, buffer_order = make_buffer_order(p, λ)

    α = p.α[buffer_order]
    Pc = p.P[buffer_order, buffer_order]
    A = make_absorbtion_matrix(p, λ)
    T = Pc[1:s,:] * A[:,2:end]    
    τ = (A'*α/sum(α))[2:end]    
    return (x)->1 - τ' * T^floor(x) * ones(s)
end


