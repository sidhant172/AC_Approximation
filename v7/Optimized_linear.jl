using PowerModels
using JuMP
using Ipopt
# using Clp
# using GLPKMathProgInterface
using Gurobi
using MAT

include("opf_mod.jl")
include("support_functions.jl")
include("find_optimal_linearization.jl")
include("find_linearization_error.jl")



################################################################################
# OPTIONS

# algorithm parameters
cnst_gen_max_iter  = 1000   # max iterations for constraint generation
# tol = 1e-4   # convergence tolerance

# operational conditions
gen_inflation = 0.2     # defining range of loading conditions
load_inflation = 0.2    # defining range of generation conditions
# v_inflation = 0.1

tol = gen_inflation*1e-2

obj_tuning = 1e2

# quantity = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"

network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("case118.m")
# network_data = PowerModels.parse_file("nesta_case57_ieee.m")
# network_data = PowerModels.parse_file("nesta_case300_ieee.m")


# line = (18,11,13)
# line = (11,5,11)


solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
# solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")
# solver_ipopt = IpoptSolver(linear_solver="ma97")

# solver_lp = GLPKSolverLP()
# solver_lp = ClpSolver()
solver_lp = GurobiSolver()



inflation_factors = Dict{String,Float64}()
inflation_factors["gen_inflation"] = gen_inflation
inflation_factors["load_inflation"] = load_inflation
# inflation_factors["v_inflation"] = v_inflation

to_approx_list = Dict{Int64,Any}()

to_approx = Dict{String,Any}()
to_approx["quantity"] = "line_real_power"
to_approx["quantity_index"] = (18,11,13)
to_approx_list[1] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "bus_voltage_magnitude"
# to_approx["quantity_index"] = 2
# to_approx_list[2] = to_approx
#
# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_reactive_power"
# to_approx["quantity_index"] = (18,11,13)
# to_approx_list[3] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_reactive_power"
# to_approx["quantity_index"] = (2,2,3)
# to_approx_list[1] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_real_power"
# to_approx["quantity_index"] = (2,2,3)
# to_approx_list[1] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_real_power"
# to_approx["quantity_index"] = (11,5,11)
# to_approx_list[1] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "bus_voltage_magnitude"
# to_approx["quantity_index"] = 5
# to_approx_list[1] = to_approx

tic()
linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()



linearation_coefficients = Dict{String,Any}()
linearation_coefficients["l0"] = linear_approximations[1]["l0"]
linearation_coefficients["l_pb"] = linear_approximations[1]["l_pb"]
linearation_coefficients["l_qb"] = linear_approximations[1]["l_qb"]
linearation_coefficients["l_v"] = linear_approximations[1]["l_v"]

(l,u) = find_linearization_error(network_data, to_approx, solver_ipopt, linearation_coefficients, inflation_factors)


########## some post checking ###########
variables = matread("ptdf_matrices.mat")
lp_jac = variables["Hac_f"][18,:][1:24]
lq_jac = variables["Hac_f"][18,:][25:48]

lp = Dict{String,Float64}()
lq = Dict{String,Float64}()

for i in network_data["active_buses"]
    lp[string(i)] = lp_jac[i]
    lq[string(i)] = lq_jac[i]
end

@show norm_error_p = sqrt( sum( (lp[string(i)] - linear_approximations[1]["l_pb"][string(i)])^2 for i in network_data["active_buses"]   )  )
@show norm_error_q = sqrt( sum( (lq[string(i)] - linear_approximations[1]["l_qb"][string(i)])^2 for i in network_data["active_buses"]   )  )
