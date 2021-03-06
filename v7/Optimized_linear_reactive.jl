using PowerModels
using JuMP
using Ipopt
# using Clp
#using GLPKMathProgInterface
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


inflation_const = 0.05
# operational conditions
gen_inflation = inflation_const    # defining range of loading conditions
load_inflation = inflation_const    # defining range of generation conditions
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

network_data_old = deepcopy(network_data)


solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
# solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")
# solver_ipopt = IpoptSolver(linear_solver="ma97")

# solver_lp = GLPKSolverLP()
# solver_lp = ClpSolver()
solver_lp = GurobiSolver(TuneOutput=0)



inflation_factors = Dict{String,Float64}()
inflation_factors["gen_inflation"] = gen_inflation
inflation_factors["load_inflation"] = load_inflation
# inflation_factors["v_inflation"] = v_inflation

to_approx_list = Dict{Int64,Any}()
# to_approx_list = Dict{String,Any}()

for (i,branch) in network_data_old["branch"]
    to_approx = Dict{String,Any}()
    to_approx["quantity"] = "line_real_power"
    to_approx["quantity_index"] = (parse(Int64,i),branch["f_bus"],branch["t_bus"])
    to_approx_list[parse(Int64,i)] = to_approx
end



tic()
linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()


# modifying linear_approximations to make it possible to write to a mat find_optimal_linearizations
num_bus = length(network_data_old["bus"])
num_branch = length(network_data["branch"])

coeff_const = zeros(num_branch)
coeff_p = zeros(num_branch,num_bus)
coeff_q = zeros(num_branch,num_bus)
approx_error = zeros(num_branch)
for (linenum,approximation) in linear_approximations
    for (i,coeff) in approximation["l_pb"]
        coeff_p[linenum,parse(Int64,i)] = coeff
    end
    for (i,coeff) in approximation["l_qb"]
        coeff_q[linenum,parse(Int64,i)] = coeff
    end
    coeff_const[linenum] = approximation["l0"]
    approx_error[linenum] = approximation["error"]
end

# write aproximations for real power
matwrite("linear_approximations_real"string(gen_inflation)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))


to_approx_list = Dict{Int64,Any}()


for (i,branch) in network_data_old["branch"]
    to_approx = Dict{String,Any}()
    to_approx["quantity"] = "line_reactive_power"
    to_approx["quantity_index"] = (parse(Int64,i),branch["f_bus"],branch["t_bus"])
    to_approx_list[parse(Int64,i)] = to_approx
end

tic()
linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()

coeff_const = zeros(num_branch)
coeff_p = zeros(num_branch,num_bus)
coeff_q = zeros(num_branch,num_bus)
approx_error = zeros(num_branch)
for (linenum,approximation) in linear_approximations
    for (i,coeff) in approximation["l_pb"]
        coeff_p[linenum,parse(Int64,i)] = coeff
    end
    for (i,coeff) in approximation["l_qb"]
        coeff_q[linenum,parse(Int64,i)] = coeff
    end
    coeff_const[linenum] = approximation["l0"]
    approx_error[linenum] = approximation["error"]
end

# write aproximations for reactive power
matwrite("results_57bus/linear_approximations_reactive"string(gen_inflation)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))
