
function make_dist( p::OverflowNetworkParameters;
                    λ = [1.5282856, 0.83785627, 1.37303848, 0.81047698, 0.83730385])
    n = p.n
    full_buffers = (1:n)[λ .> p.μ]
    s = length(full_buffers)
    s_bar = n - s
    buffer_order = vcat(full_buffers,setdiff(1:n,full_buffers))

    α = p.α[buffer_order]
    Pc = p.P[buffer_order, buffer_order]
    Qc = p.Q[buffer_order, buffer_order]
    ρ = [λ[i]/p.μ[i] for i in buffer_order]


    # B = Matrix([Diagonal([1/ρ[i] for i in 1:s]); zeros(s_bar,s)])
    
    # C = [I-Diagonal([1-1/ρ[i] for i in 1:s]) zeros(s,s_bar); zeros(s_bar,s) zeros(s_bar, s_bar)]* Qc  +  Diagonal(vcat(zeros(s),ones(s_bar))) * Pc
    
    # A = inv(Matrix(I-C))*B

    # display(A)

    α = p.α[buffer_order]
    P = p.P[buffer_order, buffer_order]
    Q = p.Q[buffer_order, buffer_order]
    γ = [p.μ[i]/λ[i] for i in buffer_order]

    C = [Q[1:s,:] .* (1 .- γ[1:s]);
         P[s+1:end, :]]
    
    pBar = 1 .- P*ones(n)
    qBar = 1 .- Q*ones(n)

    B = [((1 .- γ[1:s]) .* qBar[1:s]) Diagonal(γ[1:s]); pBar[s+1:end] zeros(s_bar,s)]
    A = inv(Matrix(I-C))*B


    T = [Matrix(I,s,s) zeros(s,s_bar)]*Pc*A[:,2:end]

    # τ = (α'*A / sum(α))
    # @show τ

    τ = make_init_dist(p)
    return (k)->1 - τ' * T^k * ones(s)
end

function make_init_dist( p::OverflowNetworkParameters;
    λ = [1.5282856, 0.83785627, 1.37303848, 0.81047698, 0.83730385])
    n = p.n
    full_buffers = (1:n)[λ .> p.μ]
    s = length(full_buffers)
    s_bar = n - s
    buffer_order = vcat(full_buffers,setdiff(1:n,full_buffers))

    α = p.α[buffer_order]
    P = p.P[buffer_order, buffer_order]
    Q = p.Q[buffer_order, buffer_order]
    γ = [p.μ[i]/λ[i] for i in buffer_order]

    C = [Q[1:s,:] .* (1 .- γ[1:s]);
         P[s+1:end, :]]
    
    pBar = 1 .- P*ones(n)
    qBar = 1 .- Q*ones(n)

    B = [((1 .- γ[1:s]) .* qBar[1:s]) Diagonal(γ[1:s]); pBar[s+1:end] zeros(s_bar,s)]
    A = inv(Matrix(I-C))*B

    println("In init...")
    display(A)

    return (A'*α/sum(α))[2:end]
end
