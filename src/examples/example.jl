
#An example 2 node network
outer_α = 2.2
example1 = OverflowNetworkParameters(   name = "2buffers",
                                        K = ones(Int,2), 
                                        α=outer_α*[0.5, 0.2],
                                        μ=[1.0, 1.7],
                                        P = [0.0 0.9; 
                                             0.2 0.0],
                                        Q = [0.0 0.1; 
                                             0.1 0.0], 
                                        α_scv = 1.0*ones(2),
                                        μ_scv = 0.5*ones(2)
                                        )
λ1 = solve_overflow_traffic_equation(example1.α, example1.μ, example1.P, example1.Q)
nfull1 = sum(λ1 - example1.μ .> 0)

#An example 5 node network
outer_α = 1.5
example2 = OverflowNetworkParameters(   name = "5buffers",
                                        K = ones(Int,5), 
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
                                        μ_scv = 0.7*ones(5)
                                        )
λ2 = solve_overflow_traffic_equation(example2.α, example2.μ, example2.P, example2.Q)
nfull2 = sum(λ2 - example2.μ .> 0)



function rand_sub_stoch_mat(L::Int , seed::Int)
    Random.seed!(seed)
    P = rand(L,L)
    P = P ./ sum(P, dims=2) #normalize rows by the sum
    P = P .* (0.2 .+ 0.7rand(L)) # multiply rows by factors in [0.2,0.9] 
    return P
end

#An example 10 node network
outer_α = 0.4
example3 = OverflowNetworkParameters(   name = "10buffers",
                                        K = ones(Int,10), 
                                        α=outer_α*ones(10),
                                        μ= ones(10),
                                        P = rand_sub_stoch_mat(10,0),
                                        Q = rand_sub_stoch_mat(10,1), 
                                        α_scv = 1.0*ones(10),
                                        μ_scv = 0.5*ones(10)
                                        )

λ3 = solve_overflow_traffic_equation(example3.α, example3.μ, example3.P, example3.Q)
nfull3 = sum(λ3 - example3.μ .> 0)

#An example 20 node network
outer_α = 0.5
example4 = OverflowNetworkParameters(   name = "20buffers",
                                        K = ones(Int,20), 
                                        α=outer_α*ones(20),
                                        μ= ones(20),
                                        P = rand_sub_stoch_mat(20,222),
                                        Q = rand_sub_stoch_mat(20,333), 
                                        α_scv = 1.0*ones(20),
                                        μ_scv = 0.5*ones(20)
                                        )

λ4 = solve_overflow_traffic_equation(example4.α, example4.μ, example4.P, example4.Q)
nfull4 = sum(λ4 - example4.μ .> 0)
