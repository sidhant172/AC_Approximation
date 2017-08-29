using PowerModels
using JuMP
using Ipopt
using Clp

include("opf_mod.jl")
include("support_functions.jl")

function find_optimal_linearizations(network_data, to_approx_list, inflation_factors, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
    ######### Defining indices and parameters ######################################
    ind_bus = [parse(Int,key) for (key,b) in network_data["bus"]]
    ind_branch = [parse(Int,key) for (key,b) in network_data["branch"]]
    ind_gen = [parse(Int,key) for (key,b) in network_data["gen"]]
    ################################################################################
    slack = [bus["index"] for (i,bus) in network_data["bus"] if bus["bus_type"] == 3][1]
    ############## Mapping generators to buses #####################################
    gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    gen_buses = unique(gen_buses)
    load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
    active_buses = union(gen_buses,load_buses)


    # remove the next 4 lines
    # gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    # gen_buses = unique(gen_buses)
    # load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
    # # active_buses = union(gen_buses,load_buses)
    # active_buses = [b["index"] for (i,b) in network_data["bus"]]

    gens_at_bus = Dict{String,Array{Int,1}}()
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data["gen"][j]["index"] for j in keys(network_data["gen"]) if network_data["gen"][j]["gen_bus"] == i]
    end
    ################################################################################
    network_data["ind_gen"] = ind_gen
    network_data["ind_bus"] = ind_bus
    network_data["ind_branch"] = ind_branch
    network_data["gen_buses"] = gen_buses
    network_data["load_buses"] = load_buses
    network_data["active_buses"] = active_buses
    network_data["gen_inflation"] = inflation_factors["gen_inflation"]
    network_data["load_inflation"] = inflation_factors["load_inflation"]
    network_data["gens_at_bus"] = gens_at_bus
    network_data["slack"] = slack

    ############# tighten generator and load limits around OPF solution based on inflation_factors ######################
    output = run_ac_opf(network_data, solver_ipopt)

    gen_inflation = inflation_factors["gen_inflation"]
    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    baseMVA = output["solution"]["baseMVA"]

    for i in ind_gen
        pg_init[i] = output["solution"]["gen"][string(i)]["pg"]/baseMVA
        qg_init[i] = output["solution"]["gen"][string(i)]["qg"]/baseMVA
        network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end

    # network_data["bus"][string(slack)]["vmax"] = (1 + inflation_factors["v_inflation"])*output["solution"]["bus"][string(slack)]["vm"]
    # network_data["bus"][string(slack)]["vmin"] = (1 - inflation_factors["v_inflation"])*output["solution"]["bus"][string(slack)]["vm"]

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

    # network_data["to_approx"] = to_approx   # add to network_data
    network_data["quantity"] = to_approx["quantity"]
    network_data["quantity_index"] = to_approx["quantity_index"]
    network_data["obj_tuning"] = obj_tuning
    active_buses = network_data["active_buses"]
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]
    slack = network_data["slack"]

################################################################################
    # m = Model(solver=ClpSolver())   # JuMP model for master problem
    # m = Model(solver=GurobiSolver())
    m = Model(solver = solver_lp)

    @variable(m, z >=0, start=0)  # maximum error
    @variable(m,-5<=l_pb[i in active_buses]<=5, start = 0)
    @variable(m,-5<=l_qb[i in active_buses]<=5, start = 0)
    @variable(m, -5<=l0<=5, start = 0) # constant term
    @variable(m, -5<=l_v<=5, start=0)

    @objective(m,Min,z)

    @constraint(m, l_pb[slack] == 0)
    @constraint(m, l_qb[slack] == 0)

    l0_val = 0
    l_v_val = 0
    l_pb_val = Dict{String,Float64}()
    l_qb_val = Dict{String,Float64}()
    for i in active_buses
        l_pb_val[string(i)] = 0.0
        l_qb_val[string(i)] = 0.0
    end


    converged_flag = 0
    #lp_err = Float64[];
    #nlp_err = Float64[];
    for iter=1:cnst_gen_max_iter   # total iterations of constraint generation
        status = solve(m)
        z_val  = getvalue(z)
        l0_val = getvalue(l0)
        l_v_val = getvalue(l_v)
        for i in active_buses
            l_pb_val[string(i)] = getvalue(l_pb[i])
            l_qb_val[string(i)] = getvalue(l_qb[i])
        end

        # add current linear coefficients to network data
        network_data["l0"] = l0_val
        network_data["l_v"] = l_v_val
        network_data["l_pb"] = l_pb_val
        network_data["l_qb"] = l_qb_val

        #push!(lp_err,abs(z_val))
        #push!(nlp_err,0)

        converged_flag = 1

    ################################################################################
        network_data["direction"] = 0   # direction of maximization
        (result, pm_model) = run_ac_opf_mod(network_data,solver)

        @show z_val - result["objective"]/obj_tuning
        #nlp_err[iter] = abs(z_val - s["objective"])

        if (z_val - result["objective"]/obj_tuning) < -tol
            converged_flag = 0

            current_sol = get_current_solution(result["solution"], pm_model, to_approx, ind_gen, ind_bus, ind_branch)
            @show current_sol["val"]

            approximation["worst_case_lower"] = current_sol

            @constraint(m, current_sol["val"] - (l0 + l_v*current_sol["vm"][slack] +
                sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) ) for i in gen_buses)
            -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
                <= z)
        end
    ################################################################################

    ################################################################################
        network_data["direction"] = 1   # direction of maximization
        (result, pm_model) = run_ac_opf_mod(network_data,solver)

        @show z_val - result["objective"]/obj_tuning
        #nlp_err[iter] = abs(z_val - s["objective"])

        if (z_val - result["objective"]/obj_tuning) < -tol
            converged_flag = 0

            current_sol = get_current_solution(result["solution"], pm_model, to_approx, ind_gen, ind_bus, ind_branch)
            @show current_sol["val"]

            approximation["worst_case_upper"] = current_sol

            @constraint(m, -current_sol["val"] + (l0 + l_v*current_sol["vm"][slack] +
            sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) ) for i in gen_buses)
        -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
            <= z)
        end
    ################################################################################

    @show iter
    @show z_val

        if converged_flag == 1
            println("Constraint generation converged!")
            break
        end
    end     # end of constraint generation loop
    ################################################################################

    approximation["l0"] = l0_val
    approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = getvalue(z)

    return approximation
end     # end of find optimal linearizaion


function find_maximum_error(network_data, to_approx, solver, obj_tuning)
    ######### Defining indices and parameters ######################################
    ind_bus = [parse(Int,key) for (key,b) in network_data["bus"]]
    ind_branch = [parse(Int,key) for (key,b) in network_data["branch"]]
    ind_gen = [parse(Int,key) for (key,b) in network_data["gen"]]
    ################################################################################
    ############## Mapping generators to buses #####################################
    gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    gen_buses = unique(gen_buses)
    load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
    active_buses = union(gen_buses,load_buses)

    # remove the next 4 lines
    # gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    # gen_buses = unique(gen_buses)
    # load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
    # # active_buses = union(gen_buses,load_buses)
    # active_buses = [b["index"] for (i,b) in network_data["bus"]]

    gens_at_bus = Dict{String,Array{Int,1}}()
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data["gen"][j]["index"] for j in keys(network_data["gen"]) if network_data["gen"][j]["gen_bus"] == i]
    end
    ################################################################################
    network_data["ind_gen"] = ind_gen
    network_data["ind_bus"] = ind_bus
    network_data["ind_branch"] = ind_branch
    network_data["gen_buses"] = gen_buses
    network_data["load_buses"] = load_buses
    network_data["active_buses"] = active_buses
    network_data["gen_inflation"] = inflation_factors["gen_inflation"]
    network_data["load_inflation"] = inflation_factors["load_inflation"]
    network_data["gens_at_bus"] = gens_at_bus

    slack = [bus["index"] for (i,bus) in network_data["bus"] if bus["bus_type"] == 3][1]
    network_data["slack"] = slack

    ############# tighten generator and load limits around OPF solution based on inflation_factors ######################
    output = run_ac_opf(network_data, solver_ipopt)

    gen_inflation = inflation_factors["gen_inflation"]
    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    baseMVA = output["solution"]["baseMVA"]

    for i in ind_gen
        pg_init[i] = output["solution"]["gen"][string(i)]["pg"]/baseMVA
        qg_init[i] = output["solution"]["gen"][string(i)]["qg"]/baseMVA
        network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end
    ################################################################################




    l0_val = 0
    l_v_val = 0
    l_pb_val = Dict{String,Float64}()
    l_qb_val = Dict{String,Float64}()
    for i in active_buses
        l_pb_val[string(i)] = 0.0
        l_qb_val[string(i)] = 0.0
    end
end
