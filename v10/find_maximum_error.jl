using PowerModels

function post_ac_opf_maxerror(data::Dict{String,Any}, model=Model(), aux_data)

    # l0_val = pm.data["l0"]
    # l_v_val = pm.data["l_v"]
    # l_pb_val = pm.data["l_pb"]
    # l_qb_val = pm.data["l_qb"]
    # obj_tuning = pm.data["obj_tuning"]

    @assert !(data["multinetwork"])
    ref = PowerModels.build_ref(data)[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    # Adding load variables
    load_inflation = aux_data["load_inflation"]
	@variable(pm.model, min((1-load_inflation)*ref[:bus][i]["pd"],(1+load_inflation)*ref[:bus][i]["pd"]) <=
        pd[i in aux_data["load_buses"]] <= max((1-load_inflation)*ref[:bus][i]["pd"],(1+load_inflation)*ref[:bus][i]["pd"]))
    @variable(pm.model, min((1-load_inflation)*ref[:bus][i]["qd"],(1+load_inflation)*ref[:bus][i]["qd"]) <=
        qd[i in aux_data["load_buses"]] <= max((1-load_inflation)*ref[:bus][i]["qd"],(1+load_inflation)*ref[:bus][i]["qd"]))

    # objective
    slack = [i for (i,bus) in ref[:ref_buses]][1]

    if aux_data["quantity"] == "line_real_power"
        @objective(model, Max, obj_tuning*(aux_data["sign"]*p[aux_data["quantity_index"] - sign*(aux_data["l0"] + aux_data["l_v"]*vm[slack] +
            sum(aux_data["l_pb"][i]*(sum(pg[j] for j in ref[:bus_gens][i]))  +  aux_data["l_qb"][i]*(sum(qg[j] for j in ref[:bus_gens][i]))   for i in aux_data["gen_buses"])
            -  sum(aux_data["l_pb"][i]*pd[i] + aux_data["l_qb"][i]*qd[i]  for i in aux_data["load_buses"])  ))
            )
        elseif aux_data["quantity"] == "line_reactive_power"
            @objective(model, Max, obj_tuning*(aux_data["sign"]*q[aux_data["quantity_index"] - sign*(aux_data["l0"] + aux_data["l_v"]*vm[slack] +
                sum(aux_data["l_pb"][i]*(sum(pg[j] for j in ref[:bus_gens][i]))  +  aux_data["l_qb"][i]*(sum(qg[j] for j in ref[:bus_gens][i]))   for i in aux_data["gen_buses"])
                -  sum(aux_data["l_pb"][i]*pd[i] + aux_data["l_qb"][i]*qd[i]  for i in aux_data["load_buses"])  ))
                )
                #done till here

        elseif aux_data["quantity"] == "bus_voltage_magnitude"
            @objective(pm.model, Max, obj_tuning*(sign*v[quantity_index] - sign*(l0_val + l_v_val*v[slack] +
            sum(l_pb_val[string(i)]*(sum(pg[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses)
            -  sum(l_pb_val[string(i)]*pd[i] + l_qb_val[string(i)]*qd[i]  for i in load_buses)  ))
            )
        else println("is not supported.")
    end
    # from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    # @objective(model, Min,
    #     sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
    #     sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    # )

    # slack bus angle
    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end

    # modifying KCL constraints to accommodate for load as variables
    for (i,bus) in ref[:bus]
        # Bus KCL
        if i in aux_data["load_buses"]
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
                sum(pg[g] for g in ref[:bus_gens][i]) - pd[i] - bus["gs"]*vm[i]^2
            )
            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
                sum(qg[g] for g in ref[:bus_gens][i]) - qd[i] + bus["bs"]*vm[i]^2
            )
        else
            # load is 0
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
                sum(pg[g] for g in ref[:bus_gens][i])  - bus["gs"]*vm[i]^2
            )
            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
                sum(qg[g] for g in ref[:bus_gens][i])  + bus["bs"]*vm[i]^2
            )
        end
    end

    # branch PF constraints
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        c = branch["br_b"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @NLconstraint(model, p_fr == g/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+c/2)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        @NLconstraint(model, p_to == g*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+c/2)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end
