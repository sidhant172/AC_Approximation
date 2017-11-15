using PowerModels
using JuMP
using Ipopt
# using Clp
#using GLPKMathProgInterface
using Gurobi
using MAT


@show filename = ARGS[1]
@show dirname = ARGS[2]
@show quantity_to_approx = ARGS[3]
@show linenum = convert(Int64,parse(Float64,ARGS[4]))
@show inflation_const = parse(Float64,ARGS[5])



include("opf_mod.jl")
include("support_functions.jl")
# include("find_optimal_linearization.jl")
include("find_optimal_linearization_gd.jl")
include("find_linearization_error.jl")



################################################################################
# OPTIONS

# algorithm parameters
cnst_gen_max_iter  = 1000   # max iterations for constraint generation
# tol = 1e-4   # convergence tolerance


# inflation_const = 0.05
# operational conditions
gen_inflation = inflation_const    # defining range of loading conditions
load_inflation = inflation_const    # defining range of generation conditions
# v_inflation = 0.1

#tol = gen_inflation*1e-2
tol = gen_inflation*1e-3

obj_tuning = 1e1

# quantity = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"

# network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("case118.m")
# network_data = PowerModels.parse_file("nesta_case57_ieee.m")
network_data = PowerModels.parse_file(filename)
# network_data = PowerModels.parse_file("nesta_case300_ieee.m")

network_data_old = deepcopy(network_data)


solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma57") # , linear_solver="ma97"
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

# for (i,branch) in network_data_old["branch"]
branch = network_data_old["branch"][string(linenum)]
to_approx = Dict{String,Any}()
to_approx["quantity"] = quantity_to_approx
to_approx["quantity_index"] = (linenum,branch["f_bus"],branch["t_bus"])
to_approx_list[linenum] = to_approx
# end



tic()
linear_approximations = find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver_ipopt, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
time = toc()



num_bus = length(network_data_old["bus"])
num_branch = length(network_data["branch"])

coeff_const = 0
coeff_p = zeros(num_bus)
coeff_q = zeros(num_bus)
approx_error = 0
approximation = linear_approximations[linenum]
for (i,coeff) in approximation["l_pb"]
    coeff_p[parse(Int64,i)] = coeff
end
for (i,coeff) in approximation["l_qb"]
    coeff_q[parse(Int64,i)] = coeff
end
coeff_const = approximation["l0"]
approx_error = approximation["error"]



# write aproximations for reactive power
if quantity_to_approx == "line_real_power"
   matwrite(string(dirname)"/linear_approximations_real"string(gen_inflation)"_line_"string(linenum)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))
elseif quantity_to_approx == "line_reactive_power"
   matwrite(string(dirname)"/linear_approximations_reactive"string(gen_inflation)"_line_"string(linenum)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))
elseif println("quantity not supported")
end
