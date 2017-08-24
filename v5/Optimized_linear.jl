using PowerModels
using JuMP
using Ipopt
using Clp
using Gurobi

include("opf_mod.jl")
include("support_functions.jl")
include("find_optimal_linearization.jl")



################################################################################
# OPTIONS

# algorithm parameters
cnst_gen_max_iter  = 1000   # max iterations for constraint generation
tol = 1e-3   # convergence tolerance

# operational conditions
gen_inflation = 0.1     # defining range of loading conditions
load_inflation = 0.1    # defining range of generation conditions


# quantity = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"

# network_data = PowerModels.parse_file("case24_ieee_rts.m")
network_data = PowerModels.parse_file("case118.m")


# line = (18,11,13)
# line = (11,5,11)


# solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")



inflation_factors = Dict{String,Float64}()
inflation_factors["gen_inflation"] = gen_inflation
inflation_factors["load_inflation"] = load_inflation

to_approx_list = Dict{Int64,Any}()
to_approx = Dict{String,Any}()
to_approx["quantity"] = "line_real_power"
to_approx["quantity_index"] = (18,11,13)
to_approx_list[1] = to_approx
to_approx = Dict{String,Any}()
to_approx["quantity"] = "bus_voltage_magnitude"
to_approx["quantity_index"] = 2
to_approx_list[2] = to_approx
to_approx = Dict{String,Any}()
to_approx["quantity"] = "line_reactive_power"
to_approx["quantity_index"] = (18,11,13)
to_approx_list[3] = to_approx
# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_real_power"
# to_approx["quantity_index"] = (11,5,11)
# to_approx_list[1] = to_approx


linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, cnst_gen_max_iter, tol)
