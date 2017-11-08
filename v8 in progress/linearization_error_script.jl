using PowerModels
using Ipopt
# using Gurobi
using GLPKMathProgInterface
using MAT

include("opf_mod.jl")
include("support_functions.jl")
include("find_linearization_error.jl")


inflation_factors = [0.05,0.1,0.2,0.3,0.4]
# filename = "case24_ieee_rts.m"

@show filename = ARGS[1]

################################################################################
# solvers
solver_ipopt = IpoptSolver(print_level=0)
# solver_lp = GurobiSolver(TuneOutput=0)
solver_lp = GLPKSolverLP()
################################################################################


# get network data and make a copy to keep the original one unmodified
network_data = PowerModels.parse_file(filename)
network_data_old = deepcopy(network_data)


num_bus = length(network_data["bus"])
num_branch = length(network_data["branch"])

gen_buses = find_gen_buses(network_data)
load_buses = find_load_buses(network_data)
active_buses = find_active_buses(network_data)
gens_at_bus = find_gens_at_bus(network_data)


pm_before = build_generic_model(network_data_old, ACPPowerModel, PowerModels.post_opf)
result_before = solve_generic_model(pm_before,solver_ipopt,solution_builder = PowerModels.get_solution)

pg_init = Dict{Int,Float64}()
qg_init = Dict{Int,Float64}()
for i in keys(network_data["gen"])
    pg_init[parse(Int64,i)] = result_before["solution"]["gen"][string(i)]["pg"]
    qg_init[parse(Int64,i)] = result_before["solution"]["gen"][string(i)]["qg"]
end

p_line_vals = Dict{String,Float64}()
q_line_vals = Dict{String,Float64}()

for (i,branch) in network_data["branch"]
    p_line_vals[i] = getvalue(pm_before.var[:p][(parse(Int64,i),branch["f_bus"],branch["t_bus"])])
    q_line_vals[i] = getvalue(pm_before.var[:q][(parse(Int64,i),branch["f_bus"],branch["t_bus"])])
end
# p_line_vals = getvalue(pm_before.var[:p][(18,11,13)])

################################################################################

################################################################################




################################################################################
################################################################################
# get Jacobian information

# read jacobian data from file
# vars = matread("ptdf_matrices.mat")
vars = matread("case"string(num_bus)"_ptdf.mat")
pp_jac = vars["Hac_f"][1:num_branch,1:num_bus]
pq_jac = vars["Hac_f"][1:num_branch,num_bus+1:end]
qp_jac = vars["Hac_f"][num_branch+1:end,1:num_bus]
qq_jac = vars["Hac_f"][num_branch+1:end,num_bus+1:end]


################################################################################
################################################################################
pb_before  = Dict{String,Float64}()
qb_before  = Dict{String,Float64}()
# finding bus injections in the base OPF case
for k in active_buses
    if k in gen_buses
        pb_before[string(k)] = sum(pg_init[j]   for j in gens_at_bus[string(k)]) - network_data_old["bus"][string(k)]["pd"]
        qb_before[string(k)] = sum(qg_init[j]   for j in gens_at_bus[string(k)]) - network_data_old["bus"][string(k)]["qd"]
    else
        pb_before[string(k)] = - network_data_old["bus"][string(k)]["pd"]
        qb_before[string(k)] = - network_data_old["bus"][string(k)]["qd"]
    end
end
################################################################################



################################################################################
# containers for error of jacobian
err_jac_real = zeros(num_branch,length(inflation_factors))
err_jac_reactive = zeros(num_branch,length(inflation_factors))

ctr = 0
for inflation_const in inflation_factors
    ctr = ctr+1
    gen_inflation = inflation_const    # defining range of loading conditions
    load_inflation = inflation_const    # defining range of generation conditions

    inflation_factors = Dict{String,Float64}()
    inflation_factors["gen_inflation"] = gen_inflation
    inflation_factors["load_inflation"] = load_inflation

    append_network_data(network_data,inflation_factors)

    define_radius_bounds(network_data, inflation_factors, pg_init, qg_init)


    for (i,branch) in network_data_old["branch"]


################################################################################
    # for real power approximation
        lp_real = Dict{String,Float64}()
        lq_real = Dict{String,Float64}()

        for j in active_buses
            lp_real[string(j)] = pp_jac[parse(Int64,i),j]
            lq_real[string(j)] = pq_jac[parse(Int64,i),j]
        end
        const_term_real = p_line_vals[i] - sum(pb_before[string(k)]*lp_real[string(k)] + qb_before[string(k)]*lq_real[string(k)] for k in active_buses)

    # for reactive power approximation
        lp_reactive = Dict{String,Float64}()
        lq_reactive = Dict{String,Float64}()

        for j in active_buses
            lp_reactive[string(j)] = qp_jac[parse(Int64,i),j]
            lq_reactive[string(j)] = qq_jac[parse(Int64,i),j]
        end
        const_term_reactive = q_line_vals[i] - sum(pb_before[string(k)]*lp_reactive[string(k)] + qb_before[string(k)]*lq_reactive[string(k)] for k in active_buses)

################################################################################

################################################################################
        # find error in real power approximation
        to_approx = Dict{String,Any}()
        to_approx["quantity"] = "line_real_power"
        to_approx["quantity_index"] = (parse(Int64,i),branch["f_bus"],branch["t_bus"])

        linearation_coefficients = Dict{String,Any}()
        linearation_coefficients["l0"] = const_term_real
        linearation_coefficients["l_pb"] = lp_real
        linearation_coefficients["l_qb"] = lq_real
        linearation_coefficients["l_v"] = 0

        @show (l,u) = find_linearization_error(network_data, to_approx, solver_ipopt, linearation_coefficients, inflation_factors)

        err_jac_real[parse(Int64,i),ctr] = max(l,u)

################################################################################

        # find error in reactive power approximation
        to_approx = Dict{String,Any}()
        to_approx["quantity"] = "line_reactive_power"
        to_approx["quantity_index"] = (parse(Int64,i),branch["f_bus"],branch["t_bus"])

        linearation_coefficients = Dict{String,Any}()
        linearation_coefficients["l0"] = const_term_reactive
        linearation_coefficients["l_pb"] = lp_reactive
        linearation_coefficients["l_qb"] = lq_reactive
        linearation_coefficients["l_v"] = 0

        @show (l,u) = find_linearization_error(network_data, to_approx, solver_ipopt, linearation_coefficients, inflation_factors)

        err_jac_reactive[parse(Int64,i),ctr] = max(l,u)


    end

    # matwrite("resultsJacobian/jacobian_error_real_"string(inflation_const)".mat",Dict("err" => err_jac_real[:,ctr]))
    # matwrite("resultsJacobian/jacobian_error_reactive_"string(inflation_const)".mat",Dict("err" => err_jac_reactive[:,ctr]))

    matwrite("resultsJacobian/case"string(num_bus)"/jacobian_error_real_"string(inflation_const)".mat",Dict("err" => err_jac_real[:,ctr]))
    matwrite("resultsJacobian/case"string(num_bus)"/jacobian_error_reactive_"string(inflation_const)".mat",Dict("err" => err_jac_reactive[:,ctr]))

end
