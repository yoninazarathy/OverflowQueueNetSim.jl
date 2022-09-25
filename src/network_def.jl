@with_kw struct OverflowNetworkParameters
    name::String = ""                               #Example name
    K::Vector{Int}                                  #The buffer sizes
    α::Vector{Float64}                              #The rates of external arrivals to queues
    μ::Vector{Float64}                              #The service rates at queues
    P::Matrix{Float64}                              #The internal routing matrix
    Q::Matrix{Float64}                              #The overflow routing matrix
    Pc::Matrix{Float64} = [P 1 .- sum(P,dims=2)]    #The completed internal routing matrix (n+1 st column is outside)
    Qc::Matrix{Float64} = [Q 1 .- sum(Q,dims=2)]    #The completed overflow routing matrix (n+1 st column is outside)
    n::Int = length(K)                              #The number of queues
    α_scv::Vector{Float64} = ones(n)                #SCVs of inter-arrival times
    μ_scv::Vector{Float64} = ones(n)                #SCVs of service times
end

mutable struct OverflowNetworkState <: State
    queues::Vector{Queue{Int}} #A vector of queues where each queue has customer ids
    params::OverflowNetworkParameters #The parameters of the queueing system
    customers::Dict{Int,Float64} #A dictionary mapping customer ids to the time of their arrival in the system
    customer_counter::Int #A unique counter of customers counting how many arrivals ever occurred and used for unique ids
    sojourn_times::Vector{Float64} #Collects Sojourn times
    OverflowNetworkState(params::OverflowNetworkParameters) = new([Queue{Int}() for _ in 1:params.n], params,Dict{Int,Float64}(),0, [])
end

queue_full(st::OverflowNetworkState, q::Int) = length(st.queues[q]) >= st.params.K[q]
queue_empty_idle(st::OverflowNetworkState, q::Int) = isempty(st.queues[q])

Base.length(net::OverflowNetworkParameters) = length(net.K)