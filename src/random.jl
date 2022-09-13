next_arrival_time(st::State, q::Int) = rand(rate_scv_gamma(st.params.α[q], st.params.α_scv[q]))

next_service_time(st::State, q::Int) = rand(rate_scv_gamma(st.params.μ[q], st.params.μ_scv[q]));

"""
Which queue to go to after service. 0 means to leave.
"""
next_queue_output_P(st::State, q::Int) = sample(1:st.params.n+1, weights(st.params.Pc[q,:])) % (st.params.n+1)

"""
Which queue to go to after overflow. 0 means to leave.
"""
next_queue_overflow_Q(st::State, q::Int) = sample(1:st.params.n+1, weights(st.params.Qc[q,:])) % (st.params.n+1)