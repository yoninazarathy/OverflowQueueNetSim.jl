"""
External arrival to the first queue
"""
struct ExternalArrivalEvent <: Event 
    q::Int #Queue to which arrival occurs
end

"""
Internal arrival either from completion of service or from overflow
"""
struct InternalArrivalEvent <: Event
    q::Int #Queue to which arrival occurs
    cid::Int #Customer id
end

"""
End of service event
"""
struct EndOfServiceAtQueueEvent <: Event
    q::Int #The index of the queue where service finished
end


function add_new_customer_in_system(time::Float64, st::OverflowNetworkState)
    st.customer_counter += 1
    st.customers[st.customer_counter] = time #label the new customer with their arrival time
    return st.customer_counter
end

function remove_customer_from_system(time::Float64, st::OverflowNetworkState, cid::Int)
    entry_time = st.customers[cid]
    sojourn_time = time - entry_time
    delete!(st.customers, cid)
    # println("Removing customer $cid, sojourn time = $sojourn_time")
    push!(st.sojourn_times, sojourn_time)
    return nothing
end

"""
Process external arrivals.
"""
function process_event(time::Float64, st::OverflowNetworkState, ev::ExternalArrivalEvent)

    cid = add_new_customer_in_system(time::Float64, st::OverflowNetworkState)

    # The next arrival from the arrival process
    new_timed_events = TimedEvent[]
    push!(new_timed_events, TimedEvent(ev, time + next_arrival_time(st, ev.q)))

    if queue_full(st, ev.q)
        destination = next_queue_overflow_Q(st, ev.q)
        if destination == 0 #if loss from system
            remove_customer_from_system(time, st, cid)
        else
            push!(new_timed_events, TimedEvent(InternalArrivalEvent(destination, cid), time))
        end
    else
        #put customer in the queue
        enqueue!(st.queues[ev.q], cid)
        
        #if this is the only job on the server engage service
        length(st.queues[ev.q]) == 1 && push!(new_timed_events,
                                                TimedEvent(EndOfServiceAtQueueEvent(ev.q), 
                                                time + next_service_time(st,ev.q)))
    end
    
    return new_timed_events
end

"""
Process internal arrivals.
"""
function process_event(time::Float64, st::State, ev::InternalArrivalEvent)

    new_timed_events = TimedEvent[]

    if queue_full(st, ev.q)
        destination = next_queue_overflow_Q(st, ev.q)
        if destination == 0 #if loss from system
            remove_customer_from_system(time, st, ev.cid)
        else
            push!(new_timed_events, TimedEvent(InternalArrivalEvent(destination, ev.cid), time))
        end
    else
        #put customer in the queue
        enqueue!(st.queues[ev.q], ev.cid)
        
        #if this is the only job on the server engage service
        length(st.queues[ev.q]) == 1 && push!(new_timed_events,
                                                TimedEvent(EndOfServiceAtQueueEvent(ev.q), 
                                                time + next_service_time(st,ev.q)))
    end
    
    return new_timed_events
end

"""
Process service completions
"""
function process_event(time::Float64, st::State, ev::EndOfServiceAtQueueEvent)

    queue = st.queues[ev.q]

    #remove the customer from the queue
    cid = dequeue!(queue)

    new_timed_events = TimedEvent[]  

    #if another customer in the queue then start a new service
    if !isempty(queue)
        push!(new_timed_events,
            TimedEvent(EndOfServiceAtQueueEvent(ev.q), time + next_service_time(st, ev.q)))   
    end

    destination = next_queue_output_P(st, ev.q)
    if destination == 0 #if departing
        remove_customer_from_system(time, st, cid)
    else #otherwise go to another queue
        push!(new_timed_events, TimedEvent(InternalArrivalEvent(destination, cid), time))
    end
    
    return new_timed_events
end


"""
This function runs the simulation and records a full trajectory.
"""
function do_experiment_traj(params::OverflowNetworkParameters; max_time, warm_up_time = 10^3)

    queues_integral = zeros(params.n)
    last_time = 0.0
    arrival_counts = zeros(Int,params.n)


    function record_queues(time::Float64, st::OverflowNetworkState) 
        (time ≥ warm_up_time) && (queues_integral += length.(st.queues)*(time-last_time)) #Use a warmup time
        last_time = time
        return nothing
    end

    function count_arrivals(time::Float64, st::OverflowNetworkState, te::TimedEvent)
        if (te.event isa InternalArrivalEvent) || (te.event isa ExternalArrivalEvent)
            arrival_counts[te.event.q] +=1
        end
    end

    state = OverflowNetworkState(params)

    do_sim( state, 
            [TimedEvent(ExternalArrivalEvent(q),0.0) for q in 1:params.n], 
            max_time = max_time, 
            general_call_back = record_queues,
            event_call_back = count_arrivals)

    return state.sojourn_times, queues_integral/max_time, arrival_counts/max_time
end