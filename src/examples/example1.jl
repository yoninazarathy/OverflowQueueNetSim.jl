#An example 5 node network
outer_α = 1.5
example1 = OverflowNetworkParameters(   K = ones(Int,5), 
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
#This is computed using https://github.com/hmjansen/overflow-algorithm
λ1 = [1.5282856, 0.83785627, 1.37303848, 0.81047698, 0.83730385]