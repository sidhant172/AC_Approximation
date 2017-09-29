include("opf_mod.jl")
include("support_functions.jl")

function find_linearization_error(network_data, to_approx, solver, linearation_coefficients,inflation_factors)

    append_network_data(network_data,inflation_factors)   # append network_data with useful data structures

    # append network_data with specific quantity being linearized
    network_data["quantity"] = to_approx["quantity"]
    network_data["quantity_index"] = to_approx["quantity_index"]
    network_data["obj_tuning"] = 100

    network_data["l0"] = linearation_coefficients["l0"]
    network_data["l_v"] = linearation_coefficients["l_v"]
    network_data["l_pb"] = linearation_coefficients["l_pb"]
    network_data["l_qb"] = linearation_coefficients["l_qb"]



    ############# tighten generator and load limits around OPF solution based on inflation_factors ######################
    # output = run_ac_opf(network_data, solver_ipopt)
    #
    # gen_inflation = network_data["gen_inflation"]
    # pg_init = Dict{Int,Float64}()
    # qg_init = Dict{Int,Float64}()
    # baseMVA = output["solution"]["baseMVA"]
    #
    # for i in network_data["ind_gen"]
    #     pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
    #     qg_init[i] = output["solution"]["gen"][string(i)]["qg"]
    #     network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
    #     network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
    #     network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    #     network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    # end
    ################################################################################


    network_data["direction"] = 0
    (result, pm) = run_ac_opf_mod(network_data,solver)
    negative_error = result["objective"]/network_data["obj_tuning"]

    network_data["direction"] = 1
    (result, pm) = run_ac_opf_mod(network_data,solver)
    positive_error = result["objective"]/network_data["obj_tuning"]

    return negative_error, positive_error
end
