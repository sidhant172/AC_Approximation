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
        val_JuMP_var = pm.var[:nw][0][:cnd][1][:p]
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "line_reactive_power"
        # val_JuMP_var = getindex(model, :q)
        val_JuMP_var = pm.var[:nw][0][:cnd][1][:q]
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "bus_voltage_magnitude"
        val = vm_curr[index]
    else println("Approximating ", quantity, " is not supported.")
    end


    @show current_sol["pg"] = pg_curr
    current_sol["qg"] = qg_curr
    current_sol["pd"] = pd_curr
    current_sol["qd"] = qd_curr
    current_sol["p"] = p_curr
    current_sol["q"] = q_curr
    @show current_sol["vm"] = vm_curr
    current_sol["va"] = va_curr
    current_sol["val"] = val

    return current_sol
end     # end of get_current_solution




# function to append network_data_scope with useful data structures so that it can be passed and accessed within all functions
function append_network_data(network_data_scope,inflation_factors)
    ######### Defining indices and parameters ######################################
    ind_bus = [parse(Int,key) for (key,b) in network_data_scope["bus"]]
    ind_branch = [parse(Int,key) for (key,b) in network_data_scope["branch"]]
    ind_gen = [parse(Int,key) for (key,b) in network_data_scope["gen"]]
    ################################################################################
    slack = [bus["index"] for (i,bus) in network_data_scope["bus"] if bus["bus_type"] == 3][1]
    ############## Mapping generators to buses #####################################
    gen_buses = [network_data_scope["gen"][i]["gen_bus"] for i in keys(network_data_scope["gen"])]
    gen_buses = unique(gen_buses)
    # load_buses = [parse(Int64,i) for i in keys(network_data_scope["bus"]) if abs(network_data_scope["bus"][i]["pd"]) + abs(network_data_scope["bus"][i]["qd"]) > 1e-2]
    load_buses = unique([load["load_bus"] for (i,load) in network_data_scope["load"]])
    active_buses = union(gen_buses,load_buses)


    num_bus = length(network_data_scope["bus"])
    bus_loads = Dict{Int64,Any}()
    for l = 1:num_bus
        bus_loads[l] = []
        for (i,load) in network_data_scope["load"]
            if load["load_bus"] == l
                push!(bus_loads[l],i)
            end
        end
    end

    for i in load_buses
        @show i
        @show network_data_scope["bus"][string(i)]["pd"] = sum(network_data_scope["load"][string(k)]["pd"] for k in bus_loads[i])
        @show network_data_scope["bus"][string(i)]["qd"] = sum(network_data_scope["load"][string(k)]["qd"] for k in bus_loads[i])
    end

    for i in setdiff(1:num_bus,load_buses)
        network_data_scope["bus"][string(i)]["pd"] = 0
        network_data_scope["bus"][string(i)]["qd"] = 0
    end


    gens_at_bus = Dict{String,Array{Int,1}}()
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data_scope["gen"][j]["index"] for j in keys(network_data_scope["gen"]) if network_data_scope["gen"][j]["gen_bus"] == i]
    end
    ################################################################################
    network_data_scope["ind_gen"] = ind_gen
    network_data_scope["ind_bus"] = ind_bus
    network_data_scope["ind_branch"] = ind_branch
    network_data_scope["gen_buses"] = gen_buses
    network_data_scope["load_buses"] = load_buses
    network_data_scope["active_buses"] = active_buses
    network_data_scope["gen_inflation"] = inflation_factors["gen_inflation"]
    network_data_scope["load_inflation"] = inflation_factors["load_inflation"]
    network_data_scope["gens_at_bus"] = gens_at_bus
    network_data_scope["slack"] = slack
end

# function to create useful structure in network data






function define_radius_bounds(network_data_scope, inflation_factors, pg_init, qg_init)
    # output = PowerModels.run_ac_opf(network_data_scope, solver)

    gen_inflation = inflation_factors["gen_inflation"]
    # pg_init = Dict{Int,Float64}()
    # qg_init = Dict{Int,Float64}()

    for i in network_data_scope["ind_gen"]
        # pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
        # qg_init[i] = output["solution"]["gen"][string(i)]["qg"]
        network_data_scope["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data_scope["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        network_data_scope["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        network_data_scope["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end
end



############## degrees_of_freedom computation #########
function find_degrees_of_freedom(network_data_scope, solver)
    gen_buses = unique(find_gen_buses(network_data_scope))
    load_buses = unique(find_load_buses(network_data_scope))
    active_buses = find_active_buses(network_data_scope)

    gen_only = setdiff(gen_buses,load_buses)
    load_only = setdiff(load_buses,gen_buses)

    result = PowerModels.run_ac_opf(network_data_scope,solver)

    gens_at_bus = Dict{String,Array{Int,1}}()
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data_scope["gen"][j]["index"] for j in keys(network_data_scope["gen"]) if network_data_scope["gen"][j]["gen_bus"] == i]
    end

    degrees_of_freedom = 0

    for i in gen_buses
        to_add_real = 0
        if abs(network_data_scope["bus"][string(i)]["pd"]) > 0
            to_add_real = 1
        end
        for j in gens_at_bus[string(i)]
            if abs(result["solution"]["gen"][string(j)]["pg"]) > 0
                to_add_real = 1
            end
        end
        degrees_of_freedom = degrees_of_freedom + to_add_real

        to_add_reactive = 0
        if abs(network_data_scope["bus"][string(i)]["qd"]) > 0
            to_add_reactive = 1
        end
        for j in gens_at_bus[string(i)]
            if abs(result["solution"]["gen"][string(j)]["qg"]) > 0
                to_add_reactive = 1
            end
        end
        degrees_of_freedom = degrees_of_freedom + to_add_reactive
    end

    for i in load_only
        if abs(network_data_scope["bus"][string(i)]["pd"]) > 0
            degrees_of_freedom = degrees_of_freedom  + 1
        end
        if abs(network_data_scope["bus"][string(i)]["qd"]) > 0
            degrees_of_freedom = degrees_of_freedom  + 1
        end
    end

    return degrees_of_freedom - 1
end


function set_warm_start(pm::GenericPowerModel ,pm_old::GenericPowerModel)
    for i in ids(pm_old,:bus)
        # voltage variables starting point
        vm = pm.var[:nw][0][:cnd][1][:vm][i]
        va = pm.var[:nw][0][:cnd][1][:va][i]
        vm_old = pm_old.var[:nw][0][:vm][i]
        va_old = pm_old.var[:nw][0][:va][i]

        # @show getvalue(vm_old)
        # @show getvalue(va_old)

        JuMP.setvalue(vm,getvalue(vm_old))
        JuMP.setvalue(va,getvalue(va_old))

        # load variables starting point
        pd = getindex(pm.model,:pd)[i]
        qd = getindex(pm.model,:qd)[i]

        pd_old = getindex(pm_old.model,:pd)[i]
        qd_old = getindex(pm_old.model,:qd)[i]

        JuMP.setvalue(pd,getvalue(pd_old))
        JuMP.setvalue(qd,getvalue(qd_old))

        # @show getvalue(pd_old)
        # @show getvalue(qd_old)
    end


    for i in ids(pm_old,:gen)
        pg = pm.var[:nw][0][:cnd][1][:pg][i]
        qg = pm.var[:nw][0][:cnd][1][:qg][i]
        pg_old = pm_old.var[:nw][0][:pg][i]
        qg_old = pm_old.var[:nw][0][:qg][i]

        JuMP.setvalue(pg,getvalue(pg_old))
        JuMP.setvalue(qg,getvalue(qg_old))
    end


    for i in ids(pm_old,:branch)
        pf = pm.var[:nw][0][:cnd][1][:p][(i,pm.data["branch"][string(i)]["f_bus"],pm.data["branch"][string(i)]["t_bus"])]
        pf_old = pm_old.var[:nw][0][:p][(i,pm.data["branch"][string(i)]["f_bus"],pm.data["branch"][string(i)]["t_bus"])]
        # @show getvalue(p_old)

        qf = pm.var[:nw][0][:cnd][1][:q][(i,pm.data["branch"][string(i)]["f_bus"],pm.data["branch"][string(i)]["t_bus"])]
        qf_old = pm_old.var[:nw][0][:q][(i,pm.data["branch"][string(i)]["f_bus"],pm.data["branch"][string(i)]["t_bus"])]
        # @show getvalue(q_old)

        JuMP.setvalue(pf,getvalue(pf_old))
        JuMP.setvalue(qf,getvalue(qf_old))




        pt = pm.var[:nw][0][:cnd][1][:p][(i,pm.data["branch"][string(i)]["t_bus"],pm.data["branch"][string(i)]["f_bus"])]
        pt_old = pm_old.var[:nw][0][:p][(i,pm.data["branch"][string(i)]["t_bus"],pm.data["branch"][string(i)]["f_bus"])]
        # @show getvalue(p_old)

        qt = pm.var[:nw][0][:cnd][1][:q][(i,pm.data["branch"][string(i)]["t_bus"],pm.data["branch"][string(i)]["f_bus"])]
        qt_old = pm_old.var[:nw][0][:q][(i,pm.data["branch"][string(i)]["t_bus"],pm.data["branch"][string(i)]["f_bus"])]
        # @show getvalue(q_old)

        JuMP.setvalue(pt,getvalue(pt_old))
        JuMP.setvalue(qt,getvalue(qt_old))
    end
end




################################################################################
################################################################################
################################################################################
# very low level support functions
function find_active_buses(network_data_scope)
    active_buses = union(find_gen_buses(network_data_scope),find_load_buses(network_data_scope))
    return active_buses
end

function find_gen_buses(network_data_scope)
    gen_buses = [network_data_scope["gen"][i]["gen_bus"] for i in keys(network_data_scope["gen"])]
    return gen_buses
end

function find_load_buses(network_data_scope)
    # load_buses = [parse(Int64,i) for i in keys(network_data_scope["bus"]) if abs(network_data_scope["bus"][i]["pd"]) + abs(network_data_scope["bus"][i]["qd"]) > 1e-2]
    load_buses = [load["load_bus"] for (i,load) in network_data_scope["load"]]

    return load_buses
end


function find_gens_at_bus(network_data_scope)
    gens_at_bus = Dict{String,Array{Int,1}}()
    gen_buses = find_gen_buses(network_data_scope)
    for i in gen_buses
        gens_at_bus[string(i)] = [network_data_scope["gen"][j]["index"] for j in keys(network_data_scope["gen"]) if network_data_scope["gen"][j]["gen_bus"] == i]
    end
    return gens_at_bus
end
