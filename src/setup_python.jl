 PyCall
function setup_python()
    py"""
    import sys
    sys.path.insert(0, "./src/src_python")
    """
    global solve_overflow_traffic_equation = pyimport("overflow_algorithm")["solve_overflow_traffic_equation"]
    return nothing
end