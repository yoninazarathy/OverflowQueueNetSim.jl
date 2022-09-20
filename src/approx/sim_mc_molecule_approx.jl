function sim_mc_for_approx_dist( p::OverflowNetworkParameters, λ::Vector{Float64})
    n = p.n
    full_buffers = (1:n)[λ .> p.μ]
    s = length(full_buffers)
    s_bar = n - s

    reward = 0.0 
    state = -sample(1:n,weights(p.α)) 
    while state !=0
        if state > 0 #inside a buffer
            if state ∈ full_buffers
                reward += p.K[state]/p.μ[state]
            end
            state = -(sample(1:(n+1), weights(p.Pc[state,:])) % (n+1))#move to an outside buffer or 0 (out) 
        else # state < 0 is an outside of a buffer
            if -state ∈ full_buffers #if buffer is an overflowing buffer
                if rand() < p.μ[-state]/λ[-state] #with given probability move in
                    state = -state #move inside the buffer means going to the positive state
                else #with complement overflow
                    state = -(sample(1:(n+1), weights(p.Qc[-state,:])) % (n+1))#move to an outside buffer or 0 (out) 
                end
            else  #buffer is not overflowing so always move in
                state = -state #move inside the buffer
            end
        end
    end

    return reward
end


mc_sim_approx_dist(p::OverflowNetworkParameters, λ::Vector{Float64}; N = 10^6) = 
                                    ecdf([sim_mc_for_approx_dist(p, λ) for _ in 1:N])
