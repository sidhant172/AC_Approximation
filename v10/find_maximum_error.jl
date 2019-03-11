using PowerModels

function post_ac_opf_maxerror(data::Dict{String,Any}, model, aux_data)
    PMs = PowerModels

    # @assert !(data["multinetwork"])
    # ref = PowerModels.build_ref(data)[:nw][0]
    ref = aux_data["ref"]

    # voltage variables
    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    # bus injection variables
    @variable(model, aux_data["pbmin"][i] <= pb[i in aux_data["active_buses"]] <= aux_data["pbmax"][i])
    @variable(model, aux_data["qbmin"][i] <= qb[i in aux_data["active_buses"]] <= aux_data["qbmax"][i])

    # branch power flow variables
    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    # dc line power flow variables
    @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    # objective
    slack = aux_data["slack"]
    if aux_data["quantity"] == "line_real_power"
        @NLobjective(model, Max, aux_data["obj_tuning"]*(aux_data["sign"]*p[aux_data["quantity_index"]] - aux_data["sign"]*(aux_data["l0"] +
            sum(aux_data["l_pb"][i]*pb[i] + aux_data["l_qb"][i]*qb[i] for i in aux_data["active_buses"]) ))
            )
    elseif aux_data["quantity"] == "line_reactive_power"
        @NLobjective(model, Max, aux_data["obj_tuning"]*(aux_data["sign"]*q[aux_data["quantity_index"]] - aux_data["sign"]*(aux_data["l0"] +
            sum(aux_data["l_pb"][i]*pb[i] + aux_data["l_qb"][i]*qb[i] for i in aux_data["active_buses"]) ))
            )
    else println("is not supported.")
    end

    # slack bus angle
    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end



    for (i,bus) in ref[:bus]
        # Bus KCL
        if i in aux_data["active_buses"]
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == pb[i] - bus["gs"]*vm[i]^2
            )

            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == qb[i] + bus["bs"]*vm[i]^2
            )
        else
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == -bus["gs"]*vm[i]^2
            )

            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == bus["bs"]*vm[i]^2
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

    var_refs = Dict{String,Any}()

    var_refs["p"] = p
    var_refs["q"] = q
    var_refs["pb"] = pb
    var_refs["qb"] = qb

    return model, var_refs

end



function post_soc_opf_maxerror(data::Dict{String,Any}, model, aux_data)
    PMs = PowerModels

    # @assert !(data["multinetwork"])
    ref = aux_data["ref"]

    # voltage squared w variables
    @variable(model, ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2, start=1.001)

    wr_min, wr_max, wi_min, wi_max = PMs.calc_voltage_product_bounds(ref[:buspairs])

    # voltage billinear w terms
    @variable(model, wr_min[bp] <= wr[bp in keys(ref[:buspairs])] <= wr_max[bp], start=1.0)
    @variable(model, wi_min[bp] <= wi[bp in keys(ref[:buspairs])] <= wi_max[bp])

    # bus injection variables
    @variable(model, aux_data["pbmin"][i] <= pb[i in aux_data["active_buses"]] <= aux_data["pbmax"][i])
    @variable(model, aux_data["qbmin"][i] <= qb[i in aux_data["active_buses"]] <= aux_data["qbmax"][i])

    # line flow variables
    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    # dc line flow variables
    @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    # objective
    slack = aux_data["slack"]
    if aux_data["quantity"] == "line_real_power"
        @objective(model, Max, aux_data["obj_tuning"]*(aux_data["sign"]*p[aux_data["quantity_index"]] - aux_data["sign"]*(aux_data["l0"] +
            sum(aux_data["l_pb"][i]*pb[i] + aux_data["l_qb"][i]*qb[i] for i in aux_data["active_buses"]) ))
            )
    elseif aux_data["quantity"] == "line_reactive_power"
        @objective(model, Max, aux_data["obj_tuning"]*(aux_data["sign"]*q[aux_data["quantity_index"]] - aux_data["sign"]*(aux_data["l0"] +
            sum(aux_data["l_pb"][i]*pb[i] + aux_data["l_qb"][i]*qb[i] for i in aux_data["active_buses"]) ))
            )
    else println("is not supported.")
    end

    # SOC constraint
    for (i,j) in keys(ref[:buspairs])
        # Voltage Product Relaxation
        @constraint(model, wr[(i,j)]^2 + wi[(i,j)]^2 <= w[i]*w[j])
    end



    for (i,bus) in ref[:bus]
        # Bus KCL
        if i in aux_data["active_buses"]
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == pb[i] - bus["gs"]*w[i]
            )

            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == qb[i] + bus["bs"]*w[i]
            )
        else
            @constraint(model,
                sum(p[a] for a in ref[:bus_arcs][i]) +
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == -bus["gs"]*w[i]
            )

            @constraint(model,
                sum(q[a] for a in ref[:bus_arcs][i]) +
                sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == bus["bs"]*w[i]
            )
        end
    end


    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]
        wr_br = wr[bp_idx]
        wi_br = wi[bp_idx]

        # Line Flow
        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        c = branch["br_b"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @constraint(model, p_fr == g/tm*w_fr + (-g*tr+b*ti)/tm*(wr_br) + (-b*tr-g*ti)/tm*(wi_br) )
        @constraint(model, q_fr == -(b+c/2)/tm*w_fr - (-b*tr-g*ti)/tm*(wr_br) + (-g*tr+b*ti)/tm*(wi_br) )

        @constraint(model, p_to == g*w_to + (-g*tr-b*ti)/tm*(wr_br) + (-b*tr+g*ti)/tm*(-wi_br) )
        @constraint(model, q_to == -(b+c/2)*w_to - (-b*tr+g*ti)/tm*(wr_br) + (-g*tr-b*ti)/tm*(-wi_br) )

        # Phase Angle Difference Limit
        @constraint(model, wi_br <= tan(branch["angmax"])*wr_br)
        @constraint(model, wi_br >= tan(branch["angmin"])*wr_br)

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

    var_refs = Dict{String,Any}()

    var_refs["p"] = p
    var_refs["q"] = q
    var_refs["pb"] = pb
    var_refs["qb"] = qb

    return model, var_refs
end
