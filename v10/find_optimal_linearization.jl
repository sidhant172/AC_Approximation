using PowerModels
using JuMP
using Ipopt
using Clp

include("find_maximum_error.jl")
include("support_functions.jl")



# function to batch process a list of approximations using outer approximation
function find_all_optimal_linearizations_outer_approximation(network_data, to_approx_list, inflation_factors, solver, solver_lp, cnst_gen_max_iter, cnst_gen_tol, obj_tuning)

    aux_data = create_aux_data(network_data,inflation_factors, solver, obj_tuning)

    linear_approximations = Dict{Int64, Any}()  # constainer for all the linear approximations
    for (i,to_approx) in to_approx_list
        aux_data["quantity"] = to_approx["quantity"]
        aux_data["quantity_index"] = to_approx["quantity_index"]
        approximation = find_optimal_linearization_outer_approximation(network_data, aux_data, solver, solver_lp, cnst_gen_max_iter, cnst_gen_tol)
        linear_approximations[i] = approximation
    end
    return linear_approximations
end

# function to batch process a list of approximations using gradient descent
function find_all_optimal_linearizations_gradient_descent(network_data, to_approx_list, inflation_factors, jacobian_filename, solver, solver_lp, max_iter, tol, obj_tuning, step_size)

    aux_data = create_aux_data(network_data,inflation_factors, solver, obj_tuning)

    linear_approximations = Dict{Int64, Any}()  # constainer for all the linear approximations
    for (i,to_approx) in to_approx_list
        aux_data["quantity"] = to_approx["quantity"]
        aux_data["quantity_index"] = to_approx["quantity_index"]
        approximation = find_optimal_linearization_gradient_descent(network_data, aux_data, jacobian_filename, solver, solver_lp, max_iter, tol, step_size)
        linear_approximations[i] = approximation
    end
    return linear_approximations
end
##########################################################################################################################################################


# function to find optimal linearization using outer approximation
function find_optimal_linearization_outer_approximation(network_data, aux_data, solver, solver_lp, cnst_gen_max_iter, cnst_gen_tol)
    approximation = Dict{String,Any}()

    ref = aux_data["ref"]
    slack = aux_data["slack"]
    gen_buses = aux_data["gen_buses"]
    load_buses = aux_data["load_buses"]
    active_buses = aux_data["active_buses"]

################################################################################
    m = Model(solver = solver_lp)

    @variable(m, z >=0, start=0)  # maximum error
    @variable(m,-1<=l_pb[i in active_buses]<=1, start = 0)
    @variable(m,-1<=l_qb[i in active_buses]<=1, start = 0)
    @variable(m, -5<=l0<=5, start = 0) # constant term
    # @variable(m, -1<=l_v<=1, start=0)

    @objective(m,Min,z)
    @constraint(m, l_pb[slack] == 0)
    # @constraint(m, l_v ==0)

    # initialize all linearization coefficients to zero
    l0_val = 0.0
    # l_v_val = 0.0
    l_pb_val = Dict{Int64,Float64}()
    l_qb_val = Dict{Int64,Float64}()
    for i in aux_data["active_buses"]
        l_pb_val[i] = 0.0
        l_qb_val[i] = 0.0
    end

    converged_flag = 0

    for iter=1:cnst_gen_max_iter   # total iterations of constraint generation
        status = solve(m)
        z_val  = getvalue(z)
        l0_val = getvalue(l0)
        # l_v_val = getvalue(l_v)
        for i in active_buses
            l_pb_val[i] = getvalue(l_pb[i])
            l_qb_val[i] = getvalue(l_qb[i])
        end

        # add to aux_data
        aux_data["l0"] = l0_val
        # aux_data["l_v"] = l_v_val
        aux_data["l_pb"] = l_pb_val
        aux_data["l_qb"] = l_qb_val

        # assume converged until proven otherwise
        converged_flag = 1
        ############################################################################
        aux_data["sign"] = 1   # direction of maximization
        model = JuMP.Model(solver = solver)
        # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
        model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
        status = solve(model)
        current_sol = get_current_solution(network_data,model,var_refs,aux_data)
        @show current_sol["quantity_val"]
        @show z_val - current_sol["objval"]

        if (z_val - current_sol["objval"]) < -cnst_gen_tol
            converged_flag = 0
            approximation["worst_case_lower"] = current_sol

            @constraint(m, current_sol["quantity_val"] - (l0 +
                sum(l_pb[i]*( sum(current_sol["pg"][j] for j in ref[:bus_gens][i]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in ref[:bus_gens][i]) ) for i in gen_buses)
            -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
                <= z)
        end
        ################################################################################
        aux_data["sign"] = -1   # direction of maximization
        model = JuMP.Model(solver = solver)
        # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
        model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
        status = solve(model)
        current_sol = get_current_solution(network_data,model,var_refs,aux_data)
        @show current_sol["quantity_val"]
        @show z_val - current_sol["objval"]

        if (z_val - current_sol["objval"]) < -cnst_gen_tol
            converged_flag = 0
            approximation["worst_case_upper"] = current_sol

            @constraint(m, -current_sol["quantity_val"] + (l0 +
                sum(l_pb[i]*( sum(current_sol["pg"][j] for j in ref[:bus_gens][i]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in ref[:bus_gens][i]) ) for i in gen_buses)
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
    # approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = getvalue(z)

    return approximation
end     # end of find optimal linearizaion


# function to find optimal linearization using gradient descent
function find_optimal_linearization_gradient_descent(network_data, aux_data, jacobian_filename, solver, solver_lp, max_iter, tol, step_size)
    approximation = Dict{String,Any}()

    # some useful quantities extracted
    active_buses = aux_data["active_buses"]
    gen_buses = aux_data["gen_buses"]
    load_buses = aux_data["load_buses"]
    num_bus = length(keys(aux_data["ref"][:bus]))
    ref = aux_data["ref"]
    slack = aux_data["slack"]

    # read the jacobian data
    vars = matread(jacobian_filename)
    pp_jac = Dict{Int,Float64}()
    pq_jac = Dict{Int,Float64}()
    qp_jac = Dict{Int,Float64}()
    qq_jac = Dict{Int,Float64}()


    if to_approx["quantity"] == "line_real_power"
        for i in active_buses
            pp_jac[i] = vars["Hac_f"][to_approx["quantity_index"][1],i]
            pq_jac[i] = vars["Hac_f"][to_approx["quantity_index"][1],num_bus+i]
        end
    elseif to_approx["quantity"] == "line_reactive_power"
        for i in active_buses
            qp_jac[i] = vars["Hac_f"][length(ind_branch)+to_approx["quantity_index"][1],i]
            qq_jac[i] = vars["Hac_f"][length(ind_branch)+to_approx["quantity_index"][1],num_bus+i]
        end
    else println("not supported")
    end



    # create linearization coefficients for the jacobian approximation
    l0_val = 0  # l0_val is always zero because of the way we solve the problem
    l_v_val = 0
    l_pb_val = Dict{Int64,Float64}()
    l_qb_val = Dict{Int64,Float64}()
    for i in active_buses
        if to_approx["quantity"] == "line_real_power"
            l_pb_val[i] = pp_jac[i]
            l_qb_val[i] = pq_jac[i]
        elseif to_approx["quantity"] == "line_reactive_power"
            l_pb_val[i] = qp_jac[i]
            l_qb_val[i] = qq_jac[i]
        end
    end

    # add linearization coefficients to aux_data
    aux_data["l0"] = l0_val
    # aux_data["l_v"] = l_v_val
    aux_data["l_pb"] = l_pb_val
    aux_data["l_qb"] = l_qb_val

    # find the error for the jacobian
    aux_data["sign"] = 1   # direction of maximization
    model = JuMP.Model(solver = solver)
    # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
    model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
    status = solve(model)
    current_sol = get_current_solution(network_data,model,var_refs,aux_data)
    pos_error = current_sol["objval"]

    aux_data["sign"] = -1   # direction of maximization
    model = JuMP.Model(solver = solver)
    # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
    model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
    status = solve(model)
    current_sol = get_current_solution(network_data,model,var_refs,aux_data)
    neg_error = current_sol["objval"]

    println("Printing error of jacobian with symmetrization ...")
    @show err = 0.5*(pos_error + neg_error)
    @show l0_by_mean = 0.5*(pos_error-neg_error)


    # gradient descent iterations
    for iter = 1:max_iter
        @show iter

        # containers for the gradient descent step initialized to zero
        step_pb = Dict{Int,Float64}()
        step_qb = Dict{Int,Float64}()
        for i in active_buses
            step_pb[i] = 0.0
            step_qb[i] = 0.0
        end


        # add linearization coefficients to aux_data
        aux_data["l0"] = l0_val
        # aux_data["l_v"] = l_v_val
        aux_data["l_pb"] = l_pb_val
        aux_data["l_qb"] = l_qb_val

        # computing gradient step
        aux_data["sign"] = 1   # direction of maximization
        model = JuMP.Model(solver = solver)
        # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
        model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
        status = solve(model)
        current_sol = get_current_solution(network_data,model,var_refs,aux_data)
        pos_error = current_sol["objval"]

        # update gradient step from positive part of the objective
        for i in gen_buses
            step_pb[i] = step_pb[i] - step_size*( sum(current_sol["pg"][j] for j in ref[:bus_gens][i]) )
            step_qb[i] = step_qb[i] - step_size*( sum(current_sol["qg"][j] for j in ref[:bus_gens][i]) )
        end
        for i in load_buses
            step_pb[i] = step_pb[i] + step_size*current_sol["pd"][i]
            step_qb[i] = step_qb[i] + step_size*current_sol["qd"][i]
        end

        aux_data["sign"] = -1   # direction of maximization
        model = JuMP.Model(solver = solver)
        # model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
        model, var_refs = post_soc_opf_maxerror(network_data, model, aux_data)
        status = solve(model)
        current_sol = get_current_solution(network_data,model,var_refs,aux_data)
        neg_error = current_sol["objval"]

        # print current error
        @show err = 0.5*(pos_error + neg_error)

        # update gradient step from negative part of the objective
        for i in gen_buses
            step_pb[i] = step_pb[i] + step_size*( sum(current_sol["pg"][j] for j in ref[:bus_gens][i]) )
            step_qb[i] = step_qb[i] + step_size*( sum(current_sol["qg"][j] for j in ref[:bus_gens][i]) )
        end
        for i in load_buses
            step_pb[i] = step_pb[i] - step_size*current_sol["pd"][i]
            step_qb[i] = step_qb[i] - step_size*current_sol["qd"][i]
        end

        # perform the gradient descent step
        for i in active_buses
            l_pb_val[i] = l_pb_val[i] - step_pb[i]
            l_qb_val[i] = l_qb_val[i] - step_qb[i]
        end
        l_pb_val[slack] = 0 # no coefficeint for active power injection at slack bus
        l0_val = 0.5*(pos_error-neg_error)

    end # end of gradient descent iterations

    approximation["l0"] = 0.5*(pos_error-neg_error)
    # approximation["l_v"] = l_v_val
    approximation["l_pb"] = l_pb_val
    approximation["l_qb"] = l_qb_val
    approximation["error"] = 0.5*(pos_error + neg_error)

    return approximation
end




function find_linearization_error(network_data, inflation_factors, to_approx, approximation, solver, obj_tuning)

    aux_data = create_aux_data(network_data,inflation_factors, solver, obj_tuning)
    aux_data["quantity"] = to_approx["quantity"]
    aux_data["quantity_index"] = to_approx["quantity_index"]
    aux_data["l0"] = approximation["l0"]
    # aux_data["l_v"] = approximation["l_v"]
    aux_data["l_pb"] = approximation["l_pb"]
    aux_data["l_qb"] = approximation["l_qb"]

    aux_data["sign"] = 1
    model = JuMP.Model(solver = solver)
    model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
    status = solve(model)
    current_sol = get_current_solution(network_data,model,var_refs,aux_data)
    pos_error = current_sol["objval"]

    aux_data["sign"] = -1
    model = JuMP.Model(solver = solver)
    model, var_refs = post_ac_opf_maxerror(network_data, model, aux_data)
    status = solve(model)
    current_sol = get_current_solution(network_data,model,var_refs,aux_data)
    neg_error = current_sol["objval"]

    return pos_error,neg_error

end
