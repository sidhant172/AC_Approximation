using PowerModels
using JuMP
using Ipopt
using Clp

include("opf_mod.jl")
include("support_functions.jl")

function find_optimal_linearizations_error_tracking(network_data, to_approx_list, inflation_factors, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)

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
        approximation = find_optimal_linearization_error_tracking(network_data, to_approx, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
        linear_approximations[i] = approximation
    end
    return linear_approximations
end     # end of find_optimal_linearizations













# find optimal linearization for a given quantity
function find_optimal_linearization_error_tracking(network_data, to_approx, solver, solver_lp, cnst_gen_max_iter, tol, obj_tuning)
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
    m = Model(solver = solver_lp)

    @variable(m, z >=0, start=0)  # maximum error
    @variable(m,-1<=l_pb[i in active_buses]<=1, start = 0)
    @variable(m,-1<=l_qb[i in active_buses]<=1, start = 0)
    @variable(m, -5<=l0<=5, start = 0) # constant term
    @variable(m, -1<=l_v<=1, start=0)

    @objective(m,Min,z)

    @constraint(m, l_pb[slack] == 0)
    # @constraint(m, l_qb[slack] == 0)
    @constraint(m, l_v ==0)

    # initialize all linearization coefficients to zero
    l0_val = 0
    l_v_val = 0
    l_pb_val = Dict{String,Float64}()
    l_qb_val = Dict{String,Float64}()
    for i in active_buses
        l_pb_val[string(i)] = 0.0
        l_qb_val[string(i)] = 0.0
    end


    converged_flag = 0
    lp_err = Float64[];
    nlp_err_pos = Float64[];
    nlp_err_neg = Float64[];
    lp_deltas = Float64[];
    nlp_deltas_pos = Float64[];
    nlp_deltas_neg = Float64[]
    num_constraint_added = Int64[];


    status = solve(m)

    l0_prev = getvalue(l0)
    l_pb_prev = Dict{String, Float64}()
    l_qb_prev = Dict{String, Float64}()
    for i in active_buses
        l_pb_prev[string(i)] = getvalue(l_pb[i])
        l_qb_prev[string(i)] = getvalue(l_qb[i])
    end


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


        deltalp = (l0_val-l0_prev)^2

        for i in active_buses
            deltalp = deltalp + (l_pb_prev[string(i)] - l_pb_val[string(i)])^2 + (l_qb_prev[string(i)] - l_qb_val[string(i)])^2
        end

        push!(lp_deltas, sqrt(deltalp))

        l0_prev = l0_val
        for i in active_buses
            l_pb_prev[string(i)] = l_pb_val[string(i)]
            l_qb_prev[string(i)] = l_qb_val[string(i)]
        end

        push!(lp_err,z_val)
        #push!(nlp_err,0)

        converged_flag = 1

        num_const = 0
    ################################################################################
        network_data["direction"] = 0   # direction of maximization
        (result, pm) = run_ac_opf_mod(network_data,solver)

        @show z_val - result["objective"]/obj_tuning
        push!(nlp_err_pos,result["objective"]/obj_tuning)
        #nlp_err[iter] = abs(z_val - s["objective"])

        if (z_val - result["objective"]/obj_tuning) < -tol
            converged_flag = 0

            num_const = num_const + 1

            current_sol = get_current_solution(result["solution"], pm, to_approx, ind_gen, ind_bus, ind_branch)
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
        (result, pm) = run_ac_opf_mod(network_data,solver)

        @show z_val - result["objective"]/obj_tuning
        push!(nlp_err_neg,result["objective"]/obj_tuning)
        #nlp_err[iter] = abs(z_val - s["objective"])

        if (z_val - result["objective"]/obj_tuning) < -tol
            converged_flag = 0

            num_const = num_const + 1

            current_sol = get_current_solution(result["solution"], pm, to_approx, ind_gen, ind_bus, ind_branch)
            @show current_sol["val"]

            approximation["worst_case_upper"] = current_sol

            @constraint(m, -current_sol["val"] + (l0 + l_v*current_sol["vm"][slack] +
            sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[string(i)]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[string(i)]) ) for i in gen_buses)
        -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
            <= z)
        end
    ################################################################################

    push!(num_constraint_added,num_const)

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
    approximation["lp_err"] = lp_err
    approximation["nlp_err_pos"] = nlp_err_pos
    approximation["nlp_err_neg"] = nlp_err_neg
    approximation["lp_deltas"] = lp_deltas
    approximation["num_constraint_added"] = num_constraint_added

    return approximation
end     # end of find optimal linearizaion
