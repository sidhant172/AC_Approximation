using PowerModels
using JuMP

PMs = PowerModels

""
function run_ac_opf_mod(file, solver; kwargs...)
    return run_opf_mod(file, ACPPowerModel, solver; kwargs...)
end

""
function run_dc_opf_mod(file, solver; kwargs...)
    return run_opf_mod(file, DCPPowerModel, solver; kwargs...)
end

""
function run_opf_mod(file, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(file, model_constructor, solver, post_opf_mod; kwargs...)
end

""


function post_opf_mod(pm::PMs.GenericPowerModel)
    PMs.variable_voltage(pm)
    PMs.variable_generation(pm)
    PMs.variable_line_flow(pm)
    PMs.variable_dcline_flow(pm)

    # Defining new variables for load real power
	@variable(pm.model, 0 <= pd[i in keys(pm.ref[:bus])] <= pm.ref[:bus][i]["pd"])

    # @objective(pm.model, Max, )
    #PMs.objective_min_fuel_cost(pm)

    PMs.constraint_voltage(pm)

    for (i,bus) in pm.ref[:ref_buses]
        PMs.constraint_theta_ref(pm, bus)
    end

    for (i,bus) in pm.ref[:bus]
        #PMs.constraint_kcl_shunt(pm, bus)
        constraint_kcl_shunt_mod(pm, bus)
    end

    for (i,branch) in pm.ref[:branch]
        PMs.constraint_ohms_yt_from(pm, branch)
        PMs.constraint_ohms_yt_to(pm, branch)

        PMs.constraint_phase_angle_difference(pm, branch)

        PMs.constraint_thermal_limit_from(pm, branch)
        PMs.constraint_thermal_limit_to(pm, branch)
    end
    for (i,dcline) in pm.ref[:dcline]
        PMs.constraint_dcline(pm, dcline)
    end
end





function constraint_kcl_shunt_mod{T <: PMs.AbstractACPForm}(pm::PMs.GenericPowerModel{T}, bus)
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
    #q_dc = getindex(pm.model, :q_dc)

    c1 = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - pd_var - gs*v^2)
    c2 = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - qd + bs*v^2)
    return Set([c1, c2])
end
