using PowerModels
using JuMP


PM = PowerModels

function run_ac_opf_mod(data, solver)
    pm = build_generic_model(data, ACPPowerModel, post_opf_mod)
    solution = solve_generic_model(pm, solver; solution_builder = PowerModels.get_solution)
    return solution, pm.model
end


function post_opf_mod(pm::PM.GenericPowerModel)
    # standard variable definitions using PowerModels inbuild functions
    PM.variable_voltage(pm)
    PM.variable_generation(pm)
    PM.variable_line_flow(pm)
    PM.variable_dcline_flow(pm)

    variable_load(pm)   # custom function to add load as variables

    add_objective(pm)   # custom function that adds objective

    add_power_factor_constraint(pm)     # custom function that adds constant power factor constraint


    # Rest of the constraints using PowerModels in-built functions
    PM.constraint_voltage(pm)
    for (i,bus) in pm.ref[:ref_buses]
        PM.constraint_theta_ref(pm, bus)
    end

    for (i,bus) in pm.ref[:bus]
        #PM.constraint_kcl_shunt(pm, bus)
        constraint_kcl_shunt_mod(pm, bus)
    end

    for (i,branch) in pm.ref[:branch]
        PM.constraint_ohms_yt_from(pm, branch)
        PM.constraint_ohms_yt_to(pm, branch)

        PM.constraint_phase_angle_difference(pm, branch)

        PM.constraint_thermal_limit_from(pm, branch)
        PM.constraint_thermal_limit_to(pm, branch)
    end
    for (i,dcline) in pm.ref[:dcline]
        PM.constraint_dcline(pm, dcline)
    end

end     # end post_opf_mod




# custom kcl constraint to accommodate loads as variables
function constraint_kcl_shunt_mod{T <: PM.AbstractACPForm}(pm::PM.GenericPowerModel{T}, bus)
    i = bus["index"]
    bus_arcs = pm.ref[:bus_arcs][i]
    bus_arcs_dc = pm.ref[:bus_arcs_dc][i]
    bus_gens = pm.ref[:bus_gens][i]

    pd = bus["pd"]
    qd = bus["qd"]
    gs = bus["gs"]
    bs = bus["bs"]

    v = getindex(pm.model, :v)[i]
    p = getindex(pm.model, :p)
    q = getindex(pm.model, :q)
    pg = getindex(pm.model, :pg)
    qg = getindex(pm.model, :qg)
    p_dc = getindex(pm.model, :p_dc)
    q_dc = getindex(pm.model, :q_dc)


    pd_var = getindex(pm.model, :pd)[i]
    qd_var = getindex(pm.model, :qd)[i]
    #q_dc = getindex(pm.model, :q_dc)

    c1 = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - pd_var - gs*v^2)
    c2 = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - qd_var + bs*v^2)
    return Set([c1, c2])
end     # end constraint_kcl_shunt_mod


# custom unction to add objective
function add_objective(pm::PM.GenericPowerModel)
    direction = pm.data["direction"]

    quantity = pm.data["quantity"]
    quantity_index = pm.data["quantity_index"]

    sign = (-1)^direction

    # access JuMP variables
    v = getindex(pm.model, :v)  # voltage magnitude variable
    t = getindex(pm.model, :t)  # voltage angle variable
    p = getindex(pm.model, :p)  # line active flow variables
    q = getindex(pm.model, :q)  # line reactive flow variables
    pg = getindex(pm.model, :pg)   # active generation variables
    qg = getindex(pm.model, :qg)   # reactive generation variables
    pd = getindex(pm.model, :pd)
    qd = getindex(pm.model, :qd)

    # get current linearization coefficients
    l0_val = pm.data["l0"]
    l_pb_val = pm.data["l_pb"]
    l_qb_val = pm.data["l_qb"]


    if quantity == "line_real_power"
        @objective(pm.model, Max, sign*p[quantity_index] - sign*(l0_val +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in gens_at_bus[i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in gens_at_bus[i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  )
            )
        elseif quantity == "line_reactive_power"
            @objective(pm.model, Max, sign*q[quantity_index] - sign*(l0_val +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in gens_at_bus[i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in gens_at_bus[i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  )
            )
        elseif quantity == "bus_voltage_magnitude"
            @objective(pm.model, Max, sign*v[quantity_index] - sign*(l0_val +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in gens_at_bus[i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in gens_at_bus[i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  )
            )
        else println("is not supported.")
    end
end

# custom function to add power factor constraint
function add_power_factor_constraint(pm::PM.GenericPowerModel)
    pd = getindex(pm.model, :pd)
    qd = getindex(pm.model, :qd)
    @constraint(pm.model, power_factor[i in keys(pm.ref[:bus])], qd[i]*pm.ref[:bus][i]["pd"] == pd[i]*pm.ref[:bus][i]["qd"])
end

# custom function to add load as variables
function variable_load(pm::GenericPowerModel)
    # Defining new variables for load real and reactive power
	@variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["pd"] <= pd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["pd"])
    @variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["qd"] <= qd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["qd"])
end
