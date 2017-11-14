using PowerModels
using Distributions

include("support_functions.jl")

function find_monte_carlo_error(network_data, to_approx_list, linearation_coefficients_list, inflation_factors, solver, num_samples)

    append_network_data(network_data,inflation_factors)   # append network_data with useful data structures

    # unpack useful data
    active_buses = network_data["active_buses"]
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]
    slack = network_data["slack"]

    ################ find region to sample #####################################
    output = run_ac_opf(network_data, solver_ipopt)

    gen_inflation = inflation_factors["gen_inflation"]
    load_inflation = inflation_factors["load_inflation"]

    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()

    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    pd_init = Dict{Int,Float64}()
    qd_init = Dict{Int,Float64}()

    pg_max = Dict{Int,Float64}()
    qg_max = Dict{Int,Float64}()
    pg_min = Dict{Int,Float64}()
    qg_min = Dict{Int,Float64}()

    pd_max = Dict{Int,Float64}()
    qd_max = Dict{Int,Float64}()
    pd_min = Dict{Int,Float64}()
    qd_min = Dict{Int,Float64}()

    for i in network_data["ind_gen"]
        pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
        qg_init[i] = output["solution"]["gen"][string(i)]["qg"]

        pg_max[i] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        pg_min[i] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        qg_max[i] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        qg_min[i] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end

    for i in network_data["ind_bus"]
        pd_max[i] = max(network_data["bus"][string(i)]["pd"]*(1+load_inflation),network_data["bus"][string(i)]["pd"]*(1-load_inflation))
        pd_min[i] = min(network_data["bus"][string(i)]["pd"]*(1+load_inflation),network_data["bus"][string(i)]["pd"]*(1-load_inflation))
        qd_max[i] = max(network_data["bus"][string(i)]["qd"]*(1+load_inflation),network_data["bus"][string(i)]["qd"]*(1-load_inflation))
        qd_min[i] = min(network_data["bus"][string(i)]["qd"]*(1+load_inflation),network_data["bus"][string(i)]["qd"]*(1-load_inflation))
    end
    ############################################################################


    # for (index,to_approx) in to_approx_list


    # create containers for monte_carlo maximums
    perr_pos_mc = Dict{Any,Float64}()
    perr_neg_mc = Dict{Any,Float64}()
    qerr_pos_mc = Dict{Any,Float64}()
    qerr_neg_mc = Dict{Any,Float64}()

    for (i,to_approx) in to_approx_list
        if to_approx["quantity"] == "line_real_power"
            perr_pos_mc[to_approx["quantity_index"]] = -100
            perr_neg_mc[to_approx["quantity_index"]] = -100
        elseif to_approx["quantity"] == "line_reactive_power"
            qerr_pos_mc[to_approx["quantity_index"]] = -100
            qerr_neg_mc[to_approx["quantity_index"]] = -100
        else
            println("Quantity not supported.")
        end
    end
    #####################################################

    for samples = 1:num_samples
        pg_samples = Dict{Int,Float64}()
        qg_samples = Dict{Int,Float64}()
        pd_samples = Dict{Int,Float64}()
        qd_samples = Dict{Int,Float64}()

        for i in network_data["ind_gen"]
            pg_samples[i] = pg_min[i] + (pg_max[i]-pg_min[i])*rand(1)[1]
            qg_samples[i] = qg_min[i] + (qg_max[i]-qg_min[i])*rand(1)[1]
            network_data["gen"][string(i)]["pg"] = pg_samples[i]
            network_data["gen"][string(i)]["qg"] = qg_samples[i]
        end
        for i in network_data["ind_bus"]
            pd_samples[i] = pd_min[i] + (pd_max[i]-pd_min[i])*rand(1)[1]
            qd_samples[i] = qd_min[i] + (qd_max[i]-qd_min[i])*rand(1)[1]
            network_data["bus"][string(i)]["pd"] = pd_samples[i]
            network_data["bus"][string(i)]["qd"] = qd_samples[i]
        end

        pm = build_generic_model(network_data, ACPPowerModel, PowerModels.post_pf)
        result = solve_generic_model(pm,solver)

        # current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)

        # result = PowerModels.run_ac_pf(network_data,solver)



        gen_buses = network_data["gen_buses"]
        load_buses = network_data["load_buses"]

        for (i,to_approx) in to_approx_list
            if to_approx["quantity"] == "line_real_power"
                p = pm.var[:nw][0][:p][to_approx["quantity_index"]]
                pval = getvalue(p)

                linearation_coefficients = linearation_coefficients_list[i]

                l0_val = linearation_coefficients["l0"]
                l_pb_val = linearation_coefficients["l_pb"]
                l_qb_val = linearation_coefficients["l_qb"]

                pval_approx = l0_val
                + sum(l_pb_val[string(i)]*(sum(pg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
                -  sum(l_pb_val[string(i)]*pd_samples[i] + l_qb_val[string(i)]*qd_samples[i]  for i in load_buses)

                perr_pos_mc[to_approx["quantity_index"]] = max(perr_pos_mc[to_approx["quantity_index"]],pval - pval_approx)
                perr_neg_mc[to_approx["quantity_index"]] = max(perr_neg_mc[to_approx["quantity_index"]],pval_approx - pval)

            elseif to_approx["quantity"] == "line_reactive_power"
                q = pm.var[:nw][0][:q][to_approx["quantity_index"]]
                qval = getvalue(q)

                linearation_coefficients = linearation_coefficients_list[i]

                l0_val = linearation_coefficients["l0"]
                l_pb_val = linearation_coefficients["l_pb"]
                l_qb_val = linearation_coefficients["l_qb"]

                qval_approx = l0_val
                + sum(l_pb_val[string(i)]*(sum(pg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
                -  sum(l_pb_val[string(i)]*pd_samples[i] + l_qb_val[string(i)]*qd_samples[i]  for i in load_buses)

                qerr_pos_mc[to_approx["quantity_index"]] = max(qerr_pos_mc[to_approx["quantity_index"]],qval - qval_approx)
                qerr_neg_mc[to_approx["quantity_index"]] = max(qerr_neg_mc[to_approx["quantity_index"]],qval_approx - qval)
            end
        end


    end
    # end of monte_carlo

    approximation_errors = Dict{Int,Any}()

    for (i,to_approx) in to_approx_list
        if to_approx["quantity"] == "line_real_power"
            errors = Dict{String,Float64}()
            errors["positive_error"] = perr_pos_mc[to_approx["quantity_index"]]
            errors["negative_error"] = perr_neg_mc[to_approx["quantity_index"]]
            errors["maximum_error"] = max(perr_pos_mc[to_approx["quantity_index"]],perr_neg_mc[to_approx["quantity_index"]])
            approximation_errors[i] = errors
        elseif to_approx["quantity"] == "line_reactive_power"
            errors = Dict{String,Float64}()
            errors["positive_error"] = qerr_pos_mc[to_approx["quantity_index"]]
            errors["negative_error"] = qerr_neg_mc[to_approx["quantity_index"]]
            errors["maximum_error"] = max(qerr_pos_mc[to_approx["quantity_index"]],qerr_neg_mc[to_approx["quantity_index"]])
            approximation_errors[i] = errors
        else println("Quantity not supported")
        end
    end

    return approximation_errors


end
