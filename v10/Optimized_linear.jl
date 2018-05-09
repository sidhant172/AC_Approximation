using PowerModels
using JuMP
using Ipopt
using Clp
using MAT
using KNITRO


include("find_maximum_error.jl")
include("support_functions.jl")
include("find_optimal_linearization.jl")


gen_inflation = inflation # defining range of loading conditions
load_inflation = inflation    # d
# load_inflation = 0.0
# v_inflation = 0.1

tol = gen_inflation*1e-3
obj_tuning = 10

# quantity = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"


network_data = PowerModels.parse_file(filename)

# solver = KnitroSolver()
solver = IpoptSolver(print_level=0, linear_solver="ma57")    #print_level=0,
# solver_lp = GLPKSolverLP()
solver_lp = ClpSolver()



inflation_factors = Dict{String,Float64}()
inflation_factors["gen_inflation"] = gen_inflation
inflation_factors["load_inflation"] = load_inflation


# create list of quantities to find linearizations for
to_approx_list = Dict{Int64,Any}()

to_approx = Dict{String,Any}()
to_approx["quantity"] = "line_real_power"
to_approx["quantity_index"] = (line_num,network_data["branch"][string(line_num)]["f_bus"],network_data["branch"][string(line_num)]["t_bus"])
to_approx_list[1] = to_approx


# testing outer approximation
cnst_gen_max_iter = 0
tic()
linear_approximations = find_all_optimal_linearizations_outer_approximation(network_data, to_approx_list, inflation_factors, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()

@show find_linearization_error(network_data,inflation_factors,to_approx_list[1],linear_approximations[1],solver,1.0)

# testing gradient descent
max_iter = 1500
step_size = 1e-3
tic()
linear_approximations = find_all_optimal_linearizations_gradient_descent(network_data, to_approx_list, inflation_factors, jacobian_filename, solver, solver_lp, max_iter, tol, obj_tuning, step_size)
time = toc()

@show find_linearization_error(network_data,inflation_factors,to_approx_list[1],linear_approximations[1],solver,1.0)
