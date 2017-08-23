using PowerModels
using JuMP


PM = PowerModels

function run_ac_opf_mod(data, solver)
    return run_opf_mod(data, ACPPowerModel, solver)
end


function run_opf_mod(data, model_constructor, solver)
    pm = build_generic_model(data, model_constructor, post_opf_mod)
    solution = solve_generic_model(pm, solver; solution_builder = PowerModels.get_solution)
    return solution, pm.model
end



function post_opf_mod(pm::PM.GenericPowerModel)
    # standard variable definitions using PowerModels inbuild functions
    PM.variable_voltage(pm)
    PM.variable_generation(pm)
    PM.variable_line_flow(pm)
    PM.variable_dcline_flow(pm)


    # Defining new variables for load real and reactive power
	@variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["pd"] <= pd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["pd"])
    @variable(pm.model, (1-load_inflation)*pm.ref[:bus][i]["qd"] <= qd[i in keys(pm.ref[:bus])] <= (1+load_inflation)*pm.ref[:bus][i]["qd"])




    # access JuMP variables
    v = getindex(pm.model, :v)  # voltage magnitude variable
    t = getindex(pm.model, :t)  # voltage angle variable
    p = getindex(pm.model, :p)  # line active flow variables
    q = getindex(pm.model, :q)  # line reactive flow variables
    pg = getindex(pm.model, :pg)   # active generation variables
    qg = getindex(pm.model, :qg)   # reactive generation variables


    # @show l0_val
    # @show l_pb_val
    # @show l_qb_val

    # Objective based on direction of maximization
    if pos == 0
        @objective(pm.model, Max, p[line] - (l0_val +
            sum(l_pb_val[i]*(sum(pg[j] for j in gens_at_bus[i]))  +  l_qb_val[i]*(sum(qg[j] for j in gens_at_bus[i]))   for i in gen_buses)
         -  sum(l_pb_val[i]*pd[i] + l_qb_val[i]*qd[i]  for i in load_buses)  )
        )
    end

    if pos == 1
        @objective(pm.model, Max, -p[line] + (l0_val +
            sum(l_pb_val[i]*(sum(pg[j] for j in gens_at_bus[i]))  +  l_qb_val[i]*(sum(qg[j] for j in gens_at_bus[i]))   for i in gen_buses)
         -  sum(l_pb_val[i]*pd[i] + l_qb_val[i]*qd[i]  for i in load_buses)  )
        )
    end


    # Additional constraints
    # constant power factor constraint
    @constraint(pm.model, power_factor[i in keys(pm.ref[:bus])], qd[i]*pm.ref[:bus][i]["pd"] == pd[i]*pm.ref[:bus][i]["qd"])
    # slack voltage constraint
    # @constraint(pm.model, v[1] == 1)    # voltage magnitude
    # @constraint(pm.model, t[1] == 0)    # voltage angle




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




# modified kcl constraint to accommodate loads as variables
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
