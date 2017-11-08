using PowerModels
using JuMP
using Ipopt
using Clp

include("opf_mod.jl")
include("support_functions.jl")

function find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)

    append_network_data(network_data,inflation_factors)   # append network_data with useful data structures

    ############# tighten generator and load limits around OPF solution based on inflation_factors ######################
    output = run_ac_opf(network_data, solver_ipopt)

    gen_inflation = inflation_factors["gen_inflation"]
    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    baseMVA = output["solution"]["baseMVA"]

    for i in network_data["ind_gen"]
        pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
        qg_init[i] = output["solution"]["gen"][string(i)]["qg"]
        network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end
    ################################################################################


    linear_approximations = Dict{Int64, Any}()  # constainer for all the linear approximations
    for (i,to_approx) in to_approx_list
        approximation = find_optimal_linearization(network_data, to_approx, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
        linear_approximations[i] = approximation
    end
    return linear_approximations
end     # end of find_optimal_linearizations













# find optimal linearization for a given quantity
function find_optimal_linearization(network_data, to_approx, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
    approximation = Dict{String,Any}()

    # append network_data with specific quantity data for this linearization
    network_data["quantity"] = to_approx["quantity"]
    network_data["quantity_index"] = to_approx["quantity_index"]
    network_data["obj_tuning"] = obj_tuning

    # extract useful data structures to use within this function
    active_buses = network_data["active_buses"]
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]
    slack = network_data["slack"]


################################################################################
# Gradient descent solver

# make this an input to the function
# vars = matread("ptdf_matrices.mat")
vars = matread("case57_ptdf.mat")
# vars = matread("case118_ptdf.mat")

pp_jac = Dict{Int,Float64}()
pq_jac = Dict{Int,Float64}()
qp_jac = Dict{Int,Float64}()
qq_jac = Dict{Int,Float64}()

@show length(ind_bus)

if to_approx["quantity"] == "line_real_power"
    for i in active_buses
        pp_jac[i] = vars["Hac_f"][to_approx["quantity_index"][1],i]
        pq_jac[i] = vars["Hac_f"][to_approx["quantity_index"][1],length(ind_bus)+i]
    end
elseif to_approx["quantity"] == "line_reactive_power"
    for i in active_buses
        qp_jac[i] = vars["Hac_f"][length(ind_branch)+to_approx["quantity_index"][1],i]
        qq_jac[i] = vars["Hac_f"][length(ind_branch)+to_approx["quantity_index"][1],length(ind_bus)+i]
    end
else println("not supported")
end

# pp_jac = vars["Hac_f"][1:length(ind_bus),1:length(ind_bus)]
# pp_jac = vars["Hac_f"][1:length(ind_bus),1:length(ind_bus)]
# pp_jac = vars["Hac_f"][1:length(ind_bus),1:length(ind_bus)]

################################################################################

l0_val = 0
l_v_val = 0
l_pb_val = Dict{String,Float64}()
l_qb_val = Dict{String,Float64}()
l_pb_val_old = Dict{String,Float64}()
l_qb_val_old = Dict{String,Float64}()
for i in active_buses
    # initialize with jacobian values for warm start
    if to_approx["quantity"] == "line_real_power"
        l_pb_val[string(i)] = pp_jac[i]
        l_qb_val[string(i)] = pq_jac[i]
        l_pb_val_old[string(i)] = pp_jac[i]
        l_qb_val_old[string(i)] = pq_jac[i]
    elseif to_approx["quantity"] == "line_reactive_power"
        l_pb_val[string(i)] = qp_jac[i]
        l_qb_val[string(i)] = qq_jac[i]
        l_pb_val_old[string(i)] = qp_jac[i]
        l_qb_val_old[string(i)] = qq_jac[i]
    end

end


step_size_const = 1e-3
val0 = 0
val1 = 0



######## solve once at the beginning to get warm start point
network_data["l0"] = l0_val
network_data["l_v"] = l_v_val
network_data["l_pb"] = l_pb_val
network_data["l_qb"] = l_qb_val


network_data["direction"] = 0   # direction of maximization
(result, pm_0_old) = run_ac_opf_mod(network_data,solver)
current_sol = get_current_solution(result["solution"], pm_0_old, to_approx, ind_gen, ind_bus, ind_branch)
val0 = result["objective"]/obj_tuning

network_data["direction"] = 1
(result, pm_1_old) = run_ac_opf_mod(network_data,solver)
current_sol = get_current_solution(result["solution"], pm_1_old, to_approx, ind_gen, ind_bus, ind_branch)
val1 = result["objective"]/obj_tuning
################################################################################

@show err = 0.5*(val0+val1)

# solver_warm = IpoptSolver(linear_solver="ma97",print_level=0,mu_init = 1e-8)
solver_warm = IpoptSolver(print_level=0,mu_init = 1e-5)
# solver_warm = IpoptSolver(print_level=0)

# solver_warm = IpoptSolver(linear_solver="ma97",mu_init = 1e-7)

step_factor = 1

warm = true
backtrack = true
ctr = 0

for iter = 1:cnst_gen_max_iter


    step_size = step_size_const/step_factor

    solver_spec = solver
    # if warm == false
    #     solver_spec = solver
    # else
    #     solver_spec = solver_warm
    # end


    @show iter

    network_data["l0"] = l0_val
    network_data["l_v"] = l_v_val
    network_data["l_pb"] = l_pb_val
    network_data["l_qb"] = l_qb_val

    step_pb = Dict{Int,Float64}()
    step_qb = Dict{Int,Float64}()
    for i in active_buses
        step_pb[i] = 0.0
        step_qb[i] = 0.0
    end

    network_data["direction"] = 0   # direction of maximization
    # (result, pm) = run_ac_opf_mod(network_data,solver)
    pm_0 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    # if warm ==  true
    #     set_warm_start(pm_0,pm_0_old)
    # end
    result = solve_generic_model(pm_0, solver_spec; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_0, to_approx, ind_gen, ind_bus, ind_branch)
    val0 = result["objective"]/obj_tuning

    for i in gen_buses
        step_pb[i] = step_pb[i] - step_size*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) )
        step_qb[i] = step_qb[i] - step_size*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) )
    end
    for i in load_buses
        step_pb[i] = step_pb[i] + step_size*current_sol["pd"][i]
        step_qb[i] = step_qb[i] + step_size*current_sol["qd"][i]
    end

    network_data["direction"] = 1
    pm_1 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    # if warm ==  true
    #     set_warm_start(pm_1,pm_1_old)
    # end
    result = solve_generic_model(pm_1, solver_spec; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    val1 = result["objective"]/obj_tuning

    for i in gen_buses
        step_pb[i] = step_pb[i] + step_size*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) )
        step_qb[i] = step_qb[i] + step_size*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) )
    end
    for i in load_buses
        step_pb[i] = step_pb[i] - step_size*current_sol["pd"][i]
        step_qb[i] = step_qb[i] - step_size*current_sol["qd"][i]
    end


    if @show 0.5*(val0 + val1) > err + 0*tol  &&  backtrack == true
        step_factor = step_factor*5
        ctr = ctr + 1
        for i in active_buses
            l_pb_val[string(i)] = l_pb_val_old[string(i)]
            l_qb_val[string(i)] = l_qb_val_old[string(i)]
        end
        warm = false
        backtrack =  false

        if ctr == 5
            break
        end
        continue # skip doing gradient descent step
    end

    backtrack = true
    warm = true

    @show err = 0.5*(val0 + val1)
    @show step_factor

    for i in active_buses
        l_pb_val_old[string(i)] = l_pb_val[string(i)]
        l_qb_val_old[string(i)] = l_qb_val[string(i)]
    end

    # Perform gradient descent step
    for i in active_buses
        l_pb_val[string(i)] = l_pb_val[string(i)] - step_pb[i]
        l_qb_val[string(i)] = l_qb_val[string(i)] - step_qb[i]
    end

    l_pb_val[string(slack)] = 0






    pm_0_old = pm_0
    pm_1_old = pm_1
end

    approximation["l0"] = l0_val
    approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = 0.5*(val0+val1)

    return approximation



# ################################################################################
#     m = Model(solver = solver_lp)
#
#     @variable(m, z >=0, start=0)  # maximum error
#     @variable(m,-1<=l_pb[i in active_buses]<=1, start = 0)
#     @variable(m,-1<=l_qb[i in active_buses]<=1, start = 0)
#     @variable(m, -5<=l0<=5, start = 0) # constant term
#     @variable(m, -1<=l_v<=1, start=0)
#
#     @objective(m,Min,z)
#
#     @constraint(m, l_pb[slack] == 0)
#     # @constraint(m, l_qb[slack] == 0)
#     @constraint(m, l_v ==0)
#
#     # initialize all linearization coefficients to zero
#     l0_val = 0
#     l_v_val = 0
#     l_pb_val = Dict{String,Float64}()
#     l_qb_val = Dict{String,Float64}()
#     for i in active_buses
#         l_pb_val[string(i)] = 0.0
#         l_qb_val[string(i)] = 0.0
#     end
#
#
#     converged_flag = 0
#     #lp_err = Float64[];
#     #nlp_err = Float64[];
#     for iter=1:cnst_gen_max_iter   # total iterations of constraint generation
#         status = solve(m)
#         z_val  = getvalue(z)
#         l0_val = getvalue(l0)
#         l_v_val = getvalue(l_v)
#         for i in active_buses
#             l_pb_val[string(i)] = getvalue(l_pb[i])
#             l_qb_val[string(i)] = getvalue(l_qb[i])
#         end
#
#         # add current linear coefficients to network data
#         network_data["l0"] = l0_val
#         network_data["l_v"] = l_v_val
#         network_data["l_pb"] = l_pb_val
#         network_data["l_qb"] = l_qb_val
#
#         #push!(lp_err,abs(z_val))
#         #push!(nlp_err,0)
#
#         converged_flag = 1
#
#     ################################################################################
#         network_data["direction"] = 0   # direction of maximization
#         (result, pm) = run_ac_opf_mod(network_data,solver)
#
#         @show z_val - result["objective"]/obj_tuning
#         #nlp_err[iter] = abs(z_val - s["objective"])
#
#         if (z_val - result["objective"]/obj_tuning) < -tol
#             converged_flag = 0
#
#             current_sol = get_current_solution(result["solution"], pm, to_approx, ind_gen, ind_bus, ind_branch)
#             @show current_sol["val"]
#
#             approximation["worst_case_lower"] = current_sol
#
#             @constraint(m, current_sol["val"] - (l0 + l_v*current_sol["vm"][slack] +
#                 sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) ) for i in gen_buses)
#             -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
#                 <= z)
#         end
#     ################################################################################
#
#     ################################################################################
#         network_data["direction"] = 1   # direction of maximization
#         (result, pm) = run_ac_opf_mod(network_data,solver)
#
#         @show z_val - result["objective"]/obj_tuning
#         #nlp_err[iter] = abs(z_val - s["objective"])
#
#         if (z_val - result["objective"]/obj_tuning) < -tol
#             converged_flag = 0
#
#             current_sol = get_current_solution(result["solution"], pm, to_approx, ind_gen, ind_bus, ind_branch)
#             @show current_sol["val"]
#
#             approximation["worst_case_upper"] = current_sol
#
#             @constraint(m, -current_sol["val"] + (l0 + l_v*current_sol["vm"][slack] +
#             sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) ) for i in gen_buses)
#         -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
#             <= z)
#         end
#     ################################################################################
#
#     @show iter
#     @show z_val
#
#         if converged_flag == 1
#             println("Constraint generation converged!")
#             break
#         end
#     end     # end of constraint generation loop
#     ################################################################################

    # approximation["l0"] = l0_val
    # approximation["l_v"] = l_v_val
    # approximation["l_pb"] = l_pb_val
    # approximation["l_qb"] = l_qb_val
    # approximation["error"] = getvalue(z)
    #
    # return approximation
end     # end of find optimal linearizaion
