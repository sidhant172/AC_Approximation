using PowerModels
using JuMP
using Ipopt
# using Clp
using GLPKMathProgInterface
# using Gurobi
using MAT


algo = 2
num_samples = 100

include("opf_mod.jl")
include("support_functions.jl")
include("monte_carlo_error.jl")
include("find_linearization_error.jl")

if algo == 0
    include("find_optimal_linearization.jl")
elseif algo == 1
    # include("find_optimal_linearization_gd_annealing.jl")
    include("find_optimal_linearization_gd.jl")
elseif algo == 2
    include("find_optimal_linearization_ipopt_gd.jl")
end

include("find_linearization_error.jl")
# include("find_optimal_linearizations_with_err.jl")


################################################################################
# OPTIONS


if algo == 0
    cnst_gen_max_iter = 1000
elseif algo == 1
    cnst_gen_max_iter = 15
elseif algo == 2
    cnst_gen_max_iter = 5
end

# algorithm parameters
# cnst_gen_max_iter  = 40   # max iterations for constraint generation
# tol = 1e-4   # convergence tolerance

# operational conditions
gen_inflation = 0.2 # defining range of loading conditions
load_inflation = 0.2    # defining range of generation conditions
# v_inflation = 0.1

tol = gen_inflation*1e-3

obj_tuning = 1e1

# quantity = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"

# network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("nesta_case14_ieee.m")
# network_data = PowerModels.parse_file("nesta_case30_as.m")
network_data = PowerModels.parse_file("case118.m")
# network_data = PowerModels.parse_file("nesta_case57_ieee.m")
# network_data = PowerModels.parse_file("nesta_case300_ieee.m")

network_data_old = deepcopy(network_data)


# line = (18,11,13)
# line = (11,5,11)


# solver_ipopt = IpoptSolver(print_level=0)#
# solver_ipopt = IpoptSolver()
# solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")
solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma57",tol=1e-12)
# solver_ipopt = IpoptSolver(linear_solver="ma97")

solver_lp = GLPKSolverLP()
# solver_lp = ClpSolver()
# solver_lp = GurobiSolver(TuneOutput=0)



inflation_factors = Dict{String,Float64}()
inflation_factors["gen_inflation"] = gen_inflation
inflation_factors["load_inflation"] = load_inflation
# inflation_factors["v_inflation"] = v_inflation

to_approx_list = Dict{Int64,Any}()


# line_num = 18
line_num = 100

to_approx = Dict{String,Any}()
to_approx["quantity"] = "line_real_power"
to_approx["quantity_index"] = (line_num,network_data["branch"][string(line_num)]["f_bus"],network_data["branch"][string(line_num)]["t_bus"])
to_approx_list[1] = to_approx

# to_approx = Dict{String,Any}()
# to_approx["quantity"] = "line_reactive_power"
# to_approx["quantity_index"] = (line_num,network_data["branch"][string(line_num)]["f_bus"],network_data["branch"][string(line_num)]["t_bus"])
# to_approx_list[2] = to_approx

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

####### remove this ########
# pm_before = build_generic_model(network_data, ACPPowerModel, PowerModels.post_opf)
# result_before = solve_generic_model(pm_before, solver_ipopt; solution_builder = PowerModels.get_solution)
# p_var = pm_before.var[:p][(18,11,13)]
# p_line_val = getvalue(p_var)
# network_data_old  = deepcopy(network_data)
############################

tic()
linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()


# linearization_coefficients_list = Dict{Int,Any}()
#
#
# for i in keys(linear_approximations)
#     linearization_coefficients = Dict{String,Any}()
#     linearization_coefficients["l0"] = linear_approximations[i]["l0"]
#     linearization_coefficients["l_pb"] = linear_approximations[i]["l_pb"]
#     linearization_coefficients["l_qb"] = linear_approximations[i]["l_qb"]
#     linearization_coefficients["l_v"] = 0
#     linearization_coefficients_list[i] = linearization_coefficients
# end



# network_data = deepcopy(network_data_old)
#
# @show find_monte_carlo_error(network_data, to_approx_list, linear_approximations, inflation_factors, solver_ipopt, num_samples)


network_data = deepcopy(network_data_old)
obj_tuning = 1
for (i,approximation) in linear_approximations
    @show find_linearization_error(network_data, to_approx, solver_ipopt, approximation,inflation_factors,obj_tuning)
end


# tic()
# linear_approximations = find_optimal_linearizations_error_tracking(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
# time = toc()

# num_bus = length(network_data_old["bus"])
# num_branch = length(network_data_old["branch"])
#
# approximation = linear_approximations[1]
# matwrite("results"string(num_bus)"/error_tracking_real_line1.mat",Dict("lp_err" => approximation["lp_err"], "nlp_err_pos" => approximation["nlp_err_pos"],
#     "nlp_err_neg" => approximation["nlp_err_neg"], "lp_deltas" => approximation["lp_deltas"]))
# approximation = linear_approximations[2]
# matwrite("results"string(num_bus)"/error_tracking_reactive_line1.mat",Dict("lp_err" => approximation["lp_err"], "nlp_err_pos" => approximation["nlp_err_pos"],
#     "nlp_err_neg" => approximation["nlp_err_neg"], "lp_deltas" => approximation["lp_deltas"]))
