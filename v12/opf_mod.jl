using PowerModels
using JuMP


PM = PowerModels

function run_ac_opf_mod(data, solver)
    pm = build_generic_model(data, ACPPowerModel, post_opf_mod)
    solution = solve_generic_model(pm, solver; solution_builder = PowerModels.get_solution)
    return solution, pm
end


function post_opf_mod(pm::PM.GenericPowerModel)
    # standard variable definitions using PowerModels inbuild functions
    PM.variable_voltage(pm, bounded = false)
    PM.variable_generation(pm, bounded = false)
    PM.variable_branch_flow(pm, bounded = false)
    PM.variable_dcline_flow(pm, bounded = false)

    variable_load(pm)   # custom function to add load as variables

    add_objective(pm)   # custom function that adds objective

    add_slack_constraint(pm)

    # add_power_factor_constraint(pm)     # custom function that adds constant power factor constraint

    data = pm.data
    # solution = output["solution"]
    # @show output = data["opf_solution_output"]
    # fix generator values
    pg = pm.var[:nw][0][:cnd][1][:pg]
    # qg = pm.var[:nw][0][:cnd][1][:qg]
    for i in ids(pm,:gen)
        gen_bus = data["gen"][string(i)]["gen_bus"]
        bus_type = data["bus"][string(gen_bus)]["bus_type"]
        if bus_type == 2
            @constraint(pm.model, pg[i] == data["opf_solution_gen"][string(i)]["pg"])
            # @constraint(pm.model, pg[i] == data["gen"][string(i)]["pmax"])
        end
    end
    # fix generator voltages
    vm = pm.var[:nw][0][:cnd][1][:vm]
    for i in ids(pm,:bus)
        if data["bus"][string(i)]["bus_type"] in [2,3]
            @show val = data["opf_solution_bus"][string(i)]["vm"]
            # val = data["bus"][string(i)]["vmax"]
            @constraint(pm.model, vm[i] == val)
        end
    end
    #
    # for (i,bus) in ref(pm, :ref_buses)
    #
    # end

    # pg = pm.var[:nw][0][:cnd][1][:pg]
    # qg = pm.var[:nw][0][:cnd][1][:qg]

    # Rest of the constraints using PowerModels in-built functions
    # PM.constraint_voltage(pm)

    for i in ids(pm,:ref_buses)
        # PM.constraint_theta_ref(pm, i)
    end

    for (i,bus) in pm.ref[:nw][0][:bus]
        #PM.constraint_kcl_shunt(pm, bus)
        constraint_kcl_shunt_mod(pm, bus)
    end

    for i in ids(pm,:branch)
        PM.constraint_ohms_yt_from(pm, i)
        PM.constraint_ohms_yt_to(pm, i)

        # PM.constraint_voltage_angle_difference(pm, i)
        #
        # PM.constraint_thermal_limit_from(pm, i)
        # PM.constraint_thermal_limit_to(pm, i)
    end
    for i in ids(pm,:dcline)
        PM.constraint_dcline(pm, i)
    end

end     # end post_opf_mod




# custom kcl constraint to accommodate loads as variables
function constraint_kcl_shunt_mod(pm::PM.GenericPowerModel, bus)

    # vm = var(pm, n, c, :vm, i)
    # p = var(pm, n, c, :p)
    # q = var(pm, n, c, :q)
    # pg = var(pm, n, c, :pg)
    # qg = var(pm, n, c, :qg)
    # p_dc = var(pm, n, c, :p_dc)
    # q_dc = var(pm, n, c, :q_dc)

    i = bus["index"]

    bus = ref(pm, 0, :bus, i)
    bus_arcs = ref(pm, 0, :bus_arcs, i)
    bus_arcs_dc = ref(pm, 0, :bus_arcs_dc, i)
    bus_gens = ref(pm, 0, :bus_gens, i)
    bus_loads = ref(pm, 0, :bus_loads, i)
    bus_shunts = ref(pm, 0, :bus_shunts, i)

    bus_pd = Dict(k => ref(pm, 0, :load, k, "pd", 1) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, 0, :load, k, "qd", 1) for k in bus_loads)

    bus_gs = Dict(k => ref(pm, 0, :shunt, k, "gs", 1) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, 0, :shunt, k, "bs", 1) for k in bus_shunts)

    # bus_arcs = pm.ref[:nw][0][:bus_arcs][i]
    # bus_arcs_dc = pm.ref[:nw][0][:bus_arcs_dc][i]
    # bus_gens = pm.ref[:nw][0][:bus_gens][i]
    # bus_gs = pm.ref[:nw][0][:bus_gs][i]

    vm = pm.var[:nw][0][:cnd][1][:vm][i]
    p = pm.var[:nw][0][:cnd][1][:p]
    q = pm.var[:nw][0][:cnd][1][:q]
    pg = pm.var[:nw][0][:cnd][1][:pg]
    qg = pm.var[:nw][0][:cnd][1][:qg]
    p_dc = pm.var[:nw][0][:cnd][1][:p_dc]
    q_dc = pm.var[:nw][0][:cnd][1][:q_dc]

    pd_var = getindex(pm.model, :pd)[i]
    qd_var = getindex(pm.model, :qd)[i]

    @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - pd_var  - sum(gs for gs in values(bus_gs))*vm^2)
    @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - qd_var  + sum(bs for bs in values(bus_bs))*vm^2)

    ###################################################################################

    # i = bus["index"]
    # bus_arcs = pm.ref[:nw][0][:bus_arcs][i]
    # bus_arcs_dc = pm.ref[:nw][0][:bus_arcs_dc][i]
    # bus_gens = pm.ref[:nw][0][:bus_gens][i]
    #
    # pd = bus["pd"]
    # qd = bus["qd"]
    # gs = bus["gs"]
    # bs = bus["bs"]

    # v = getindex(pm.model, :v)[i]
    # p = getindex(pm.model, :p)
    # q = getindex(pm.model, :q)
    # pg = getindex(pm.model, :pg)
    # qg = getindex(pm.model, :qg)
    # p_dc = getindex(pm.model, :p_dc)
    # q_dc = getindex(pm.model, :q_dc)

    # v = pm.var[:nw][0][:cnd][1][:vm][i]
    # p = pm.var[:nw][0][:cnd][1][:p]
    # q = pm.var[:nw][0][:cnd][1][:q]
    # pg = pm.var[:nw][0][:cnd][1][:pg]
    # qg = pm.var[:nw][0][:cnd][1][:qg]
    # p_dc = pm.var[:nw][0][:cnd][1][:p_dc]
    # q_dc = pm.var[:nw][0][:cnd][1][:q_dc]



    # pd_var = getindex(pm.model, :pd)[i]
    # qd_var = getindex(pm.model, :qd)[i]

    # pd_var = pm.var[:pd][i]
    # qd_var = pm.var[:qd][i]


    # c1 = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - pd_var - gs*v^2)
    # c2 = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - qd_var + bs*v^2)
    # return Set([c1, c2])
end     # end constraint_kcl_shunt_mod


# custom unction to add objective
function add_objective(pm::PM.GenericPowerModel)
    direction = pm.data["direction"]

    quantity = pm.data["quantity"]
    quantity_index = pm.data["quantity_index"]

    gen_buses = pm.data["gen_buses"]
    load_buses = pm.data["load_buses"]
    active_buses = pm.data["active_buses"]
    # gens_at_bus = pm.data["gens_at_bus"]
    # bus_gens = pm.ref[:bus_gens][i]

    sign = (-1)^direction

    # @show pm.var[:nw][0][:cnd][1][:cnd][1]
    v = pm.var[:nw][0][:cnd][1][:vm]
    t = pm.var[:nw][0][:cnd][1][:va]
    p = pm.var[:nw][0][:cnd][1][:p]
    q = pm.var[:nw][0][:cnd][1][:q]
    pg = pm.var[:nw][0][:cnd][1][:pg]
    qg = pm.var[:nw][0][:cnd][1][:qg]
    pd = getindex(pm.model, :pd)
    qd = getindex(pm.model, :qd)
    # pd = pm.var[:pd]
    # qd = pm.var[:qd]

    # get current linearization coefficients
    l0_val = pm.data["l0"]
    l_v_val = pm.data["l_v"]
    l_pb_val = pm.data["l_pb"]
    l_qb_val = pm.data["l_qb"]
    obj_tuning = pm.data["obj_tuning"]

    slack = pm.data["slack"]


    if quantity == "line_real_power"
        @objective(pm.model, Max, obj_tuning*(sign*p[quantity_index] - sign*(l0_val + l_v_val*v[slack] +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  ))
            )
        elseif quantity == "line_reactive_power"
            @objective(pm.model, Max, obj_tuning*(sign*q[quantity_index] - sign*(l0_val + l_v_val*v[slack] +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  ))
            )
        elseif quantity == "bus_voltage_magnitude"
            @objective(pm.model, Max, obj_tuning*(sign*v[quantity_index] - sign*(l0_val + l_v_val*v[slack] +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  ))
            )
        else println("is not supported.")
    end
end

# custom function to add power factor constraint
function add_power_factor_constraint(pm::PM.GenericPowerModel)
    pd = getindex(pm.model, :pd)
    qd = getindex(pm.model, :qd)
    # pd = pm.var[:pd]
    # qd = pm.var[:qd]
    @constraint(pm.model, power_factor[i in keys(pm.ref[:bus])], qd[i]*pm.ref[:bus][i]["pd"] == pd[i]*pm.ref[:bus][i]["qd"])
    # @constraint(pm.model, power_factor_low[i in keys(pm.ref[:bus])], 0.95*pd[i]*pm.ref[:bus][i]["qd"] <= qd[i]*pm.ref[:bus][i]["pd"])
    # @constraint(pm.model, power_factor_high[i in keys(pm.ref[:bus])], qd[i]*pm.ref[:bus][i]["pd"] <= 1.05*pd[i]*pm.ref[:bus][i]["qd"])
end

# custom function to add load as variables
function variable_load(pm::GenericPowerModel)
    # Defining new variables for load real and reactive power
    @show load_inflation = pm.data["load_inflation"]
    load_buses = pm.data["load_buses"]


    # for i in ids(pm,:bus)
    #     @show i
    #     @show pm.data["bus"][string(i)]["pd"], pm.data["bus"][string(i)]["qd"]
    # end
	# @variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["pd"] <= pd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["pd"])
    # @variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["qd"] <= qd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["qd"])
    # @variable(pm.model, min((1-load_inflation)*pm.data["bus"][string(i)]["pd"],(1+load_inflation)*pm.data["bus"][string(i)]["pd"]) <=
    #     pd[i in keys(ids(pm,:bus))] <= max((1-load_inflation)*pm.data["bus"][string(i)]["pd"],(1+load_inflation)*pm.data["bus"][string(i)]["pd"]))
    # @variable(pm.model, min((1-load_inflation)*pm.data["bus"][string(i)]["qd"],(1+load_inflation)*pm.data["bus"][string(i)]["qd"]) <=
    #     qd[i in keys(ids(pm,:bus))] <= max((1-load_inflation)*pm.ref[:bus][i]["qd"],(1+load_inflation)*ids(pm,:bus)["qd"]))

    @variable(pm.model, min((1-load_inflation)*pm.data["bus"][string(i)]["pd"],(1+load_inflation)*pm.data["bus"][string(i)]["pd"]) <=
        pd[i in ids(pm,:bus)] <= max((1-load_inflation)*pm.data["bus"][string(i)]["pd"],(1+load_inflation)*pm.data["bus"][string(i)]["pd"]))
    @variable(pm.model, min((1-load_inflation)*pm.data["bus"][string(i)]["qd"],(1+load_inflation)*pm.data["bus"][string(i)]["qd"]) <=
        qd[i in ids(pm,:bus)] <= max((1-load_inflation)*pm.data["bus"][string(i)]["qd"],(1+load_inflation)*pm.data["bus"][string(i)]["qd"]))
end

function add_slack_constraint(pm::GenericPowerModel)
    # t = getindex(pm.model, :t)
    t = pm.var[:nw][0][:cnd][1][:va]
    @constraint(pm.model, t[pm.data["slack"]] == 0)
end
