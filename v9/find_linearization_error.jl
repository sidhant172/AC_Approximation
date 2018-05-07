include("opf_mod.jl")
include("support_functions.jl")

function find_linearization_error(network_data, to_approx, solver, linearation_coefficients,inflation_factors,obj_tuning)

    append_network_data(network_data,inflation_factors)   # append network_data with useful data structures

    ############# tighten generator and load limits around OPF solution based on inflation_factors ######################
    solver_spec = IpoptSolver(print_level=0, linear_solver="ma57",tol=1e-12)
    network_data_copy = deepcopy(network_data)
    output = run_ac_opf(network_data_copy, solver_spec)

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
    ##############################################################################################################



    # append network_data with specific quantity being linearized
    network_data["quantity"] = to_approx["quantity"]
    network_data["quantity_index"] = to_approx["quantity_index"]
    network_data["obj_tuning"] = obj_tuning

    network_data["l0"] = linearation_coefficients["l0"]
    network_data["l_v"] = linearation_coefficients["l_v"]
    network_data["l_pb"] = linearation_coefficients["l_pb"]
    network_data["l_qb"] = linearation_coefficients["l_qb"]



    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]

    @show val_nominal = network_data["l0"] + sum(network_data["l_pb"][string(i)]*(sum(pg_init[j] for j in network_data["gens_at_bus"][string(i)]))  +  network_data["l_qb"][string(i)]*(sum(qg_init[j] for j in network_data["gens_at_bus"][string(i)]))   for i in gen_buses) - sum(network_data["l_pb"][string(i)]*network_data["bus"][string(i)]["pd"] + network_data["l_qb"][string(i)]*network_data["bus"][string(i)]["qd"]  for i in load_buses)

    network_data_copy = deepcopy(network_data)
    network_data_copy["direction"] = 0
    solver_spec = IpoptSolver(print_level=0, linear_solver="ma57",tol=1e-12)
    (result, pm) = run_ac_opf_mod(network_data_copy,solver_spec)
    negative_error = result["objective"]/network_data["obj_tuning"]

    network_data_copy = deepcopy(network_data)
    network_data_copy["direction"] = 1
    solver_spec = IpoptSolver(print_level=0, linear_solver="ma57",tol=1e-12)
    (result, pm) = run_ac_opf_mod(network_data_copy,solver_spec)
    positive_error = result["objective"]/network_data["obj_tuning"]

    return negative_error, positive_error
end
