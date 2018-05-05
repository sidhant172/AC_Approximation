using PowerModels

# create aux_data with additional useful structures
function create_aux_data(network_data,inflation_factors, solver, obj_tuning)
    @assert !(network_data["multinetwork"])
    # PowerModels reference objects containing all network data information
    ref = PowerModels.build_ref(network_data)[:nw][0]

    ############## Mapping generators to buses #####################################
    gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    gen_buses = unique(gen_buses)
    load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-6]
    active_buses = union(gen_buses,load_buses)

    # find load limits using load inflation_factors
    pdmin = Dict{Int64,Float64}()
    pdmax = Dict{Int64,Float64}()
    qdmin = Dict{Int64,Float64}()
    qdmax = Dict{Int64,Float64}()

    load_inflation = inflation_factors["load_inflation"]
    for i in load_buses
        pdmin[i] = min((1-load_inflation)*ref[:bus][i]["pd"],(1+load_inflation)*ref[:bus][i]["pd"])
        pdmax[i] = max((1-load_inflation)*ref[:bus][i]["pd"],(1+load_inflation)*ref[:bus][i]["pd"])
        qdmin[i] = min((1-load_inflation)*ref[:bus][i]["qd"],(1+load_inflation)*ref[:bus][i]["qd"])
        qdmax[i] = max((1-load_inflation)*ref[:bus][i]["qd"],(1+load_inflation)*ref[:bus][i]["qd"])
    end

    # find gen limits using gen inflation_factors
    pgmin = Dict{Int64,Float64}()
    pgmax = Dict{Int64,Float64}()
    qgmin = Dict{Int64,Float64}()
    qgmax = Dict{Int64,Float64}()

    result = PowerModels.run_ac_opf(network_data, solver)

    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    gen_inflation = inflation_factors["gen_inflation"]
    for i in keys(ref[:gen])
        pg_init[i] = result["solution"]["gen"][string(i)]["pg"]
        qg_init[i] = result["solution"]["gen"][string(i)]["qg"]
        pgmax[i] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        pgmin[i] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        qgmax[i] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        qgmin[i] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end
    #######

    slack = [i for (i,bus) in ref[:ref_buses]][1]

    aux_data = Dict{String,Any}()
    aux_data["gen_buses"] = gen_buses
    aux_data["load_buses"] = load_buses
    aux_data["active_buses"] = active_buses
    aux_data["pgmin"] = pgmin
    aux_data["pgmax"] = pgmax
    aux_data["qgmin"] = qgmin
    aux_data["qgmax"] = qgmax
    aux_data["pdmin"] = pdmin
    aux_data["pdmax"] = pdmax
    aux_data["qdmin"] = qdmin
    aux_data["qdmax"] = qdmax
    aux_data["slack"] = slack
    aux_data["obj_tuning"] = obj_tuning
    aux_data["ref"] = ref

    return aux_data
end




# get current get_current_solution
function get_current_solution(network_data, model, var_refs, aux_data)
    current_sol = Dict{String,Any}()

    current_sol["objval"] = getobjectivevalue(model)/aux_data["obj_tuning"]
    ref = aux_data["ref"]

    pg_curr = Dict{Int,Float64}()
    qg_curr = Dict{Int,Float64}()
    # vm_curr = Dict{Int,Float64}()
    # va_curr = Dict{Int,Float64}()
    p_curr = Dict{Any,Float64}()
    q_curr = Dict{Any,Float64}()
    pd_curr = Dict{Int,Float64}()
    qd_curr = Dict{Int,Float64}()

    for i in keys(ref[:gen])
        pg_curr[i] = getvalue(var_refs["pg"][i])
        qg_curr[i] = getvalue(var_refs["qg"][i])
    end
    for i in aux_data["load_buses"]
        pd_curr[i] = getvalue(var_refs["pd"][i])
        qd_curr[i] = getvalue(var_refs["qd"][i])
    end
    # for i in keys(ref[:bus])
    #     vm_curr[i] = getvalue(var_refs["vm"][i])
    #     va_curr[i] = getvalue(var_refs["va"][i])
    # end

    for (l,i,j) in ref[:arcs]
        p_curr[(l,i,j)] = getvalue(var_refs["p"][(l,i,j)])
        q_curr[(l,i,j)] = getvalue(var_refs["q"][(l,i,j)])
    end

    quantity_val = 0
    if aux_data["quantity"] == "line_real_power"
        quantity_val = getvalue(var_refs["p"][aux_data["quantity_index"]])
    elseif quantity == "line_reactive_power"
        quantity_val = getvalue(var_refs["q"][aux_data["quantity_index"]])
    else println("Approximating ", aux_data["quantity"], " is not supported.")
    end

    current_sol["pg"] = pg_curr
    current_sol["qg"] = qg_curr
    current_sol["pd"] = pd_curr
    current_sol["qd"] = qd_curr
    current_sol["p"] = p_curr
    current_sol["q"] = q_curr
    # current_sol["vm"] = vm_curr
    # current_sol["va"] = va_curr
    current_sol["quantity_val"] = quantity_val

    return current_sol
end
