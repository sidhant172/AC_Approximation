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



    function f(x...)
        l_pb = Dict{String,Float64}()
        l_qb = Dict{String,Float64}()

        ctr = 0
        for i in active_buses
            ctr = ctr+1
            l_pb[string(i)] = x[ctr]
        end
        for i in active_buses
            ctr = ctr+1
            l_qb[string(i)] = x[ctr]
        end

        return feval(l_pb,l_qb,network_data,solver)
    end

    function gradf(g,x...)
        l_pb = Dict{String,Float64}()
        l_qb = Dict{String,Float64}()

        ctr = 0
        for i in active_buses
            ctr = ctr+1
            l_pb[string(i)] = x[ctr]
        end
        for i in active_buses
            ctr = ctr+1
            l_qb[string(i)] = x[ctr]
        end

        gradfeval(g,l_pb,l_qb,network_data,solver)
    end



    m = Model(solver=IpoptSolver(linear_solver="ma57",tol=1e-5,hessian_approximation="limited-memory",max_iter=cnst_gen_max_iter,warm_start_init_point="yes",limited_memory_max_history=10,limited_memory_max_skipping=1,dual_inf_tol=5.0,watchdog_shortened_iter_trigger=20))
    # m = Model(solver=IpoptSolver(linear_solver="ma57",tol=1e-5,hessian_approximation="limited-memory",acceptable_obj_change_tol=1e-3,acceptable_iter=3,acceptable_dual_inf_tol=1e3))

    JuMP.register(m,:f,2*length(active_buses),f,gradf)
    @variable(m,x[1:2*length(active_buses)])
    ctr = 0
    for i in active_buses
        ctr = ctr+1
        setvalue(x[ctr],l_pb_val[string(i)])
    end
    for i in active_buses
        ctr = ctr+1
        setvalue(x[ctr],l_qb_val[string(i)])
    end

    JuMP.setNLobjective(m, :Min, Expr(:call, :f, x...))

    status = solve(m)


    ctr = 0
    for i in active_buses
        ctr = ctr+1
        l_pb_val[string(i)] = getvalue(x[ctr])
    end
    for i in active_buses
        ctr = ctr+1
        l_qb_val[string(i)] = getvalue(x[ctr])
    end



    approximation["l0"] = l0_eval(l_pb_val,l_qb_val,network_data,solver)
    approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = 0.5*(val0+val1)

    return approximation

end     # end of find optimal linearizaion


function l0_eval(l_pb,l_qb,network_data,solver)
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]

    network_data["l0"] = 0
    network_data["l_v"] = 0
    network_data["l_pb"] = l_pb
    network_data["l_qb"] = l_qb

    network_data["direction"] = 0
    pm_0 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_0, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_0, to_approx, ind_gen, ind_bus, ind_branch)
    val0 = result["objective"]/obj_tuning

    network_data["direction"] = 1
    pm_1 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_1, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    val1 = result["objective"]/obj_tuning

    return 0.5*(val0-val1)
end


function feval(l_pb,l_qb,network_data,solver)
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    active_buses = network_data["active_buses"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]

    network_data["l0"] = 0
    network_data["l_v"] = 0
    network_data["l_pb"] = l_pb
    network_data["l_qb"] = l_qb

    network_data["direction"] = 0
    pm_0 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_0, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_0, to_approx, ind_gen, ind_bus, ind_branch)
    val0 = result["objective"]/obj_tuning

    network_data["direction"] = 1
    pm_1 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_1, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    val1 = result["objective"]/obj_tuning

    return 0.5*(val0+val1)
end

function gradfeval(g,l_pb,l_qb,network_data,solver)

    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    active_buses = network_data["active_buses"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]

    g_pb = Dict{Int,Float64}()
    g_qb = Dict{Int,Float64}()
    for i in active_buses
        g_pb[i] = 0.0
        g_qb[i] = 0.0
    end

    network_data["l0"] = 0
    network_data["l_v"] = 0
    network_data["l_pb"] = l_pb
    network_data["l_qb"] = l_qb


    network_data["direction"] = 0
    pm_0 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_0, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_0, to_approx, ind_gen, ind_bus, ind_branch)
    val0 = result["objective"]/obj_tuning

    for i in gen_buses
        g_pb[i] = g_pb[i] - sum(current_sol["pg"][j] for j in gens_at_bus[string(i)])
        g_qb[i] = g_qb[i] - sum(current_sol["qg"][j] for j in gens_at_bus[string(i)])
    end
    for i in load_buses
        g_pb[i] = g_pb[i] + current_sol["pd"][i]
        g_qb[i] = g_qb[i] + current_sol["qd"][i]
    end

    network_data["direction"] = 1
    pm_1 = build_generic_model(network_data, ACPPowerModel, post_opf_mod)
    result = solve_generic_model(pm_1, solver; solution_builder = PowerModels.get_solution)
    current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)
    val1 = result["objective"]/obj_tuning

    for i in gen_buses
        g_pb[i] = g_pb[i] + sum(current_sol["pg"][j] for j in gens_at_bus[string(i)])
        g_qb[i] = g_qb[i] + sum(current_sol["qg"][j] for j in gens_at_bus[string(i)])
    end
    for i in load_buses
        g_pb[i] = g_pb[i] - current_sol["pd"][i]
        g_qb[i] = g_qb[i] - current_sol["qd"][i]
    end

    g_pb[network_data["slack"]] = 0



    ctr = 0
    for i in network_data["active_buses"]
        ctr = ctr+1
        g[ctr] = g_pb[i]
    end
    for i in network_data["active_buses"]
        ctr = ctr+1
        g[ctr] = g_qb[i]
    end

end
