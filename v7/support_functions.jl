# SUPPORT FUNCTIONS
function get_current_solution(solution::Dict{String,Any}, pm::GenericPowerModel, to_approx, ind_gen, ind_bus, ind_branch)
    baseMVA = solution["baseMVA"]

    current_sol = Dict{String,Any}()

    pg_curr = Dict{Int,Float64}()
    qg_curr = Dict{Int,Float64}()
    vm_curr = Dict{Int,Float64}()
    va_curr = Dict{Int,Float64}()
    p_curr = Dict{Array{Int,1},Float64}()
    q_curr = Dict{Array{Int,1},Float64}()
    pd_curr = Dict{Int,Float64}()
    qd_curr = Dict{Int,Float64}()



    pd = getindex(pm.model, :pd)
    qd = getindex(pm.model, :qd)
    # pd = pm.var[:pd]
    # qd = pm.var[:qd]


    for i in ind_gen
        # pg_curr[i] = solution["gen"][string(i)]["pg"]/baseMVA
        # qg_curr[i] = solution["gen"][string(i)]["qg"]/baseMVA
        pg_curr[i] = solution["gen"][string(i)]["pg"]
        qg_curr[i] = solution["gen"][string(i)]["qg"]
    end

    for i in ind_bus
        vm_curr[i] = solution["bus"][string(i)]["vm"]
        pd_curr[i] = getvalue(pd[i])
        qd_curr[i] = getvalue(qd[i])
    end


    quantity = to_approx["quantity"]
    index = to_approx["quantity_index"]
    val = 0
    if quantity == "line_real_power"
        # val_JuMP_var = getindex(model, :p)
        val_JuMP_var = pm.var[:p]
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "line_reactive_power"
        # val_JuMP_var = getindex(model, :q)
        val_JuMP_var = pm.var[:q]
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "bus_voltage_magnitude"
        val = vm_curr[index]
    else println("Approximating ", quantity, " is not supported.")
    end


    current_sol["pg"] = pg_curr
    current_sol["qg"] = qg_curr
    current_sol["pd"] = pd_curr
    current_sol["qd"] = qd_curr
    current_sol["p"] = p_curr
    current_sol["q"] = q_curr
    current_sol["vm"] = vm_curr
    current_sol["va"] = va_curr
    current_sol["val"] = val

    return current_sol
end     # end of get_current_solution




# function to append network_data with useful data structures so that it can be passed and accessed within all functions
function append_network_data(network_data,inflation_factors)
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
end

# function to create useful structure in network data






function define_radius_bounds(network_data, inflation_factors, pg_init, qg_init)
    # output = PowerModels.run_ac_opf(network_data, solver)

    gen_inflation = inflation_factors["gen_inflation"]
    # pg_init = Dict{Int,Float64}()
    # qg_init = Dict{Int,Float64}()

    for i in network_data["ind_gen"]
        # pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
        # qg_init[i] = output["solution"]["gen"][string(i)]["qg"]
        network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end
end





################################################################################
################################################################################
################################################################################
# very low level support functions
function find_active_buses(network_data)
    active_buses = union(find_gen_buses(network_data),find_load_buses(network_data))
    return active_buses
end

function find_gen_buses(network_data)
    gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
    return gen_buses
end

function find_load_buses(network_data)
    load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
    return load_buses
end


function find_gens_at_bus(network_data)
    gens_at_bus = Dict{String,Array{Int,1}}()
    gen_buses = find_gen_buses(network_data)
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data["gen"][j]["index"] for j in keys(network_data["gen"]) if network_data["gen"][j]["gen_bus"] == i]
    end
    return gens_at_bus
end
