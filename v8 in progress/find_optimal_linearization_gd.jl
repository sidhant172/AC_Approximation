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

vars = matread("case"string(length(ind_bus))"_ptdf.mat")
# make this an input to the function
# vars = matread("ptdf_matrices.mat")
# vars = matread("case14_ptdf.mat")
# vars = matread("case57_ptdf.mat")
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
for trial = 1:5
    for i in ind_bus
        vm_var = pm_0_old.var[:nw][0][:vm][i]
        setvalue(vm_var, 1 + 0.02*(2*rand()-1))
    end
    result = solve_generic_model(pm_0_old,solver)
    current_sol_temp = get_current_solution(result["solution"], pm_0_old, to_approx, ind_gen, ind_bus, ind_branch)
    if val0 < result["objective"]/obj_tuning
        val0 = result["objective"]/obj_tuning
        current_sol = deepcopy(current_sol_temp)
    end
end

network_data["direction"] = 1
(result, pm_1_old) = run_ac_opf_mod(network_data,solver)
current_sol = get_current_solution(result["solution"], pm_1_old, to_approx, ind_gen, ind_bus, ind_branch)
val1 = result["objective"]/obj_tuning
for trial = 1:5
    for i in ind_bus
        vm_var = pm_1_old.var[:nw][0][:vm][i]
        setvalue(vm_var, 1 + 0.02*(2*rand()-1))
    end
    result = solve_generic_model(pm_1_old,solver)
    current_sol_temp = get_current_solution(result["solution"], pm_1_old, to_approx, ind_gen, ind_bus, ind_branch)
    if val1 < result["objective"]/obj_tuning
        val1 = result["objective"]/obj_tuning
        current_sol = deepcopy(current_sol_temp)
    end
end

################################################################################

println("printing error of jacobian with symmetrization")
@show err = 0.5*(val0+val1)
@show l0_by_mean = 0.5*(val0-val1)

# solver_warm = IpoptSolver(linear_solver="ma97",print_level=0,mu_init = 1e-8)
solver_warm = IpoptSolver(print_level=0,mu_init = 1e-5)
# solver_warm = IpoptSolver(print_level=0)

# solver_warm = IpoptSolver(linear_solver="ma97",mu_init = 1e-7)

step_factor = 1

warm = true
backtrack = true
ctr = 0
l0_by_mean = 0

############# gradient descent iterations #############
for iter = 1:cnst_gen_max_iter


    step_size = step_size_const/step_factor

    solver_spec = solver
    # if warm == false
    #     solver_spec = solver
    # else
    #     solver_spec = solver_warm
    # end


    @show iter

    network_data["l0"] = 0
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
    pm_0 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_0,solver_spec; solution_builder = PowerModels.get_solution)
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
    result = solve_generic_model(pm_1,solver_spec; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    val1 = result["objective"]/obj_tuning

    # result = solve_generic_model(pm_1, solver_spec; solution_builder = PowerModels.get_solution)
    # current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    # val1 = result["objective"]/obj_tuning

    for i in gen_buses
        step_pb[i] = step_pb[i] + step_size*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) )
        step_qb[i] = step_qb[i] + step_size*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) )
    end
    for i in load_buses
        step_pb[i] = step_pb[i] - step_size*current_sol["pd"][i]
        step_qb[i] = step_qb[i] - step_size*current_sol["qd"][i]
    end

    @show step_norm = sqrt(sum(step_pb[i]^2 + step_qb[i]^2 for i in active_buses))
    if @show step_norm < 1e-6
        break
    end

    # if @show 0.5*(val0 + val1) > err   &&  backtrack == true
    if @show 0.5*(val0 + val1) > err + 1e-4
        @show step_factor = step_factor*2
        # ctr = ctr + 1
        for i in active_buses
            # l_pb_val[string(i)] = l_pb_val_old[string(i)] + 1e-4*(2*rand()-1)
            # l_qb_val[string(i)] = l_qb_val_old[string(i)] + 1e-4*(2*rand()-1)
            l_pb_val[string(i)] = l_pb_val_old[string(i)]
            l_qb_val[string(i)] = l_qb_val_old[string(i)]
        end
        warm = false
        # backtrack =  false


        # if @show sqrt(sum(step_pb[i]^2 + step_qb[i]^2 for i in active_buses)) < 1e-6
        #     break
        # end
        # if ctr == 7
        #     break
        # end
        continue # skip doing gradient descent step
    end

    backtrack = true
    warm = true

    @show err = 0.5*(val0 + val1)
    @show l0_by_mean = 0.5*(val0-val1)
    @show step_factor

    for i in active_buses
        l_pb_val_old[string(i)] = l_pb_val[string(i)]
        l_qb_val_old[string(i)] = l_qb_val[string(i)]
    end

    # Perform gradient descent step
    # @show sqrt(sum(step_pb[i]^2 + step_qb[i]^2 for i in active_buses))
    for i in active_buses
        l_pb_val[string(i)] = l_pb_val[string(i)] - step_pb[i]
        l_qb_val[string(i)] = l_qb_val[string(i)] - step_qb[i]
    end

    l_pb_val[string(slack)] = 0


    pm_0_old = pm_0
    pm_1_old = pm_1

end


@show     approximation["l0"] = l0_by_mean
    approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = err

    return approximation



end     # end of find optimal linearizaion
