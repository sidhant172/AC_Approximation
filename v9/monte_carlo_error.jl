using PowerModels
using Distributions

include("support_functions.jl")

function post_pf_mod(pm::GenericPowerModel)
    # PowerModels.variable_voltage(pm, bounded = false)
    # PowerModels.variable_generation(pm, bounded = false)
    # PowerModels.variable_branch_flow(pm, bounded = false)
    # PowerModels.variable_dcline_flow(pm, bounded = false)

    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm)
    PowerModels.variable_dcline_flow(pm)

    PowerModels.constraint_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        PowerModels.constraint_theta_ref(pm, i)
        # constraint_voltage_magnitude_setpoint(pm, i)
    end

    for (i,bus) in ref(pm, :bus)

        PowerModels.constraint_kcl_shunt(pm, i)

        if length(ref(pm, :bus_gens, i)) > 0
            # this assumes inactive generators are filered out of bus_gens
            # @assert bus["bus_type"] == 2

            # PowerModels.constraint_voltage_magnitude_setpoint(pm, i)
            if !(i in ids(pm,:ref_buses))
                for j in ref(pm, :bus_gens, i)
                    # PowerModels.constraint_active_gen_setpoint(pm, j)
                    pg_var = pm.var[:nw][0][:qg][j]
                    @constraint(pm.model,pg_var == pm.data["gen"][string(j)]["pg"])

                    qg_var = pm.var[:nw][0][:qg][j]
                    @constraint(pm.model,qg_var == pm.data["gen"][string(j)]["qg"])
                end
            else
                for j in ref(pm, :bus_gens, i)
                    qg_var = pm.var[:nw][0][:qg][j]
                    @constraint(pm.model,qg_var == pm.data["gen"][string(j)]["qg"])
                end
            end
        end


        # PV Bus Constraints
        # if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
        #     # this assumes inactive generators are filtered out of bus_gens
        #     @assert bus["bus_type"] == 2
        #
        #     PowerModels.constraint_voltage_magnitude_setpoint(pm, i)
        #     for j in ref(pm, :bus_gens, i)
        #         PowerModels.constraint_active_gen_setpoint(pm, j)
        #     end
        # end
    end

    for i in ids(pm, :branch)
        PowerModels.constraint_ohms_yt_from(pm, i)
        PowerModels.constraint_ohms_yt_to(pm, i)

        PowerModels.constraint_voltage_angle_difference(pm, i)
        PowerModels.constraint_thermal_limit_from(pm, i)
        PowerModels.constraint_thermal_limit_to(pm, i)
    end

    for (i,dcline) in ref(pm, :dcline)
        #constraint_dcline(pm, i) not needed, active power flow fully defined by dc line setpoints
        PowerModels.constraint_active_dcline_setpoint(pm, i)

        f_bus = ref(pm, :bus)[dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            PowerModels.constraint_voltage_magnitude_setpoint(pm, f_bus["index"])
        end

        t_bus = ref(pm, :bus)[dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            PowerModels.constraint_voltage_magnitude_setpoint(pm, t_bus["index"])
        end
    end
end



function find_monte_carlo_error(network_data, to_approx_list, linearization_coefficients_list, inflation_factors, solver, num_samples)

    pm = build_generic_model(network_data,ACPPowerModel,PowerModels.post_opf)
    output = solve_generic_model(pm,solver)
    # output = PowerModels.run_ac_opf(network_data, solver)

    qvar = pm.var[:nw][0][:q][(30,17,18)]
    @show getvalue(qvar)

    append_network_data(network_data,inflation_factors)   # append network_data with useful data structures

    # unpack useful data
    active_buses = network_data["active_buses"]
    gen_buses = network_data["gen_buses"]
    load_buses = network_data["load_buses"]
    gens_at_bus = network_data["gens_at_bus"]
    ind_gen = network_data["ind_gen"]
    ind_bus = network_data["ind_bus"]
    ind_branch = network_data["ind_branch"]
    slack = network_data["slack"]

    ################ find region to sample #####################################

    gen_inflation = inflation_factors["gen_inflation"]
    load_inflation = inflation_factors["load_inflation"]

    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()

    pg_init = Dict{Int,Float64}()
    qg_init = Dict{Int,Float64}()
    pd_init = Dict{Int,Float64}()
    qd_init = Dict{Int,Float64}()

    pg_max = Dict{Int,Float64}()
    qg_max = Dict{Int,Float64}()
    pg_min = Dict{Int,Float64}()
    qg_min = Dict{Int,Float64}()

    pd_max = Dict{Int,Float64}()
    qd_max = Dict{Int,Float64}()
    pd_min = Dict{Int,Float64}()
    qd_min = Dict{Int,Float64}()

    for i in network_data["ind_gen"]
        pg_init[i] = output["solution"]["gen"][string(i)]["pg"]
        qg_init[i] = output["solution"]["gen"][string(i)]["qg"]

        pg_max[i] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        pg_min[i] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
        qg_max[i] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
        qg_min[i] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    end

    for i in network_data["ind_bus"]
        pd_max[i] = max(network_data["bus"][string(i)]["pd"]*(1+load_inflation),network_data["bus"][string(i)]["pd"]*(1-load_inflation))
        pd_min[i] = min(network_data["bus"][string(i)]["pd"]*(1+load_inflation),network_data["bus"][string(i)]["pd"]*(1-load_inflation))
        qd_max[i] = max(network_data["bus"][string(i)]["qd"]*(1+load_inflation),network_data["bus"][string(i)]["qd"]*(1-load_inflation))
        qd_min[i] = min(network_data["bus"][string(i)]["qd"]*(1+load_inflation),network_data["bus"][string(i)]["qd"]*(1-load_inflation))
    end
    ############################################################################
    @show pg_min
    @show pg_max
    @show qg_min
    @show qg_max
    @show pd_min
    @show pd_max
    @show qd_min
    @show qd_max

    # for (index,to_approx) in to_approx_list


    # create containers for monte_carlo maximums
    perr_pos_mc = Dict{Any,Float64}()
    perr_neg_mc = Dict{Any,Float64}()
    qerr_pos_mc = Dict{Any,Float64}()
    qerr_neg_mc = Dict{Any,Float64}()

    for (i,to_approx) in to_approx_list
        if to_approx["quantity"] == "line_real_power"
            perr_pos_mc[to_approx["quantity_index"]] = -100
            perr_neg_mc[to_approx["quantity_index"]] = -100
        elseif to_approx["quantity"] == "line_reactive_power"
            qerr_pos_mc[to_approx["quantity_index"]] = -100
            qerr_neg_mc[to_approx["quantity_index"]] = -100
        else
            println("Quantity not supported.")
        end
    end
    #####################################################




    for samples = 1:num_samples
        pg_samples = Dict{Int,Float64}()
        qg_samples = Dict{Int,Float64}()
        pd_samples = Dict{Int,Float64}()
        qd_samples = Dict{Int,Float64}()

        for i in network_data["ind_gen"]
            pg_samples[i] = pg_min[i] + (pg_max[i]-pg_min[i])*rand()
            qg_samples[i] = qg_min[i] + (qg_max[i]-qg_min[i])*rand()

            network_data["gen"][string(i)]["pg"] = pg_min[i] + (pg_max[i]-pg_min[i])*rand()
            network_data["gen"][string(i)]["qg"] = qg_min[i] + (qg_max[i]-qg_min[i])*rand()
        end
        for i in network_data["ind_bus"]
            pd_samples[i] = pd_min[i] + (pd_max[i]-pd_min[i])*rand()
            qd_samples[i] = qd_min[i] + (qd_max[i]-qd_min[i])*rand()

            network_data["bus"][string(i)]["pd"] = pd_min[i] + (pd_max[i]-pd_min[i])*rand()
            network_data["bus"][string(i)]["qd"] = qd_min[i] + (qd_max[i]-qd_min[i])*rand()
        end

        # pm = build_generic_model(network_data, ACPPowerModel, post_pf_mod)
        # result = solve_generic_model(pm,solver)
        # if result["status"] == :infeasible
        #     println("Infeasible")
        #     continue
        # end


        pm = build_generic_model(network_data, ACPPowerModel, PowerModels.post_opf)

        for (i,bus) in ref(pm, :bus)
            if length(ref(pm, :bus_gens, i)) > 0
                if !(i in ids(pm,:ref_buses))
                    for j in ref(pm, :bus_gens, i)
                        # PowerModels.constraint_active_gen_setpoint(pm, j)
                        pg_var = pm.var[:nw][0][:pg][j]
                        @constraint(pm.model,pg_var == pm.data["gen"][string(j)]["pg"])

                        qg_var = pm.var[:nw][0][:qg][j]
                        @constraint(pm.model,qg_var == pm.data["gen"][string(j)]["qg"])
                    end
                else
                    for j in ref(pm, :bus_gens, i)
                        qg_var = pm.var[:nw][0][:qg][j]
                        @constraint(pm.model,qg_var == pm.data["gen"][string(j)]["qg"])
                    end
                end
            end
        end

        result = solve_generic_model(pm,solver)
        if result["status"] == :infeasible
            println("Infeasible")
            continue
        end


        # current_sol = get_current_solution(result["solution"], pm_1, to_approx, ind_gen, ind_bus, ind_branch)

        # result = PowerModels.run_ac_pf(network_data,solver)



        gen_buses = network_data["gen_buses"]
        load_buses = network_data["load_buses"]

        for (approx_num,to_approx) in to_approx_list
            if to_approx["quantity"] == "line_real_power"
                pvar = pm.var[:nw][0][:p][to_approx["quantity_index"]]
                @show pval = getvalue(pvar)

                linearization_coefficients = linearization_coefficients_list[approx_num]

                l0_val = linearization_coefficients["l0"]
                l_pb_val = linearization_coefficients["l_pb"]
                l_qb_val = linearization_coefficients["l_qb"]

                pval_approx = l0_val + sum(l_pb_val[string(i)]*(sum(pg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses) - sum(l_pb_val[string(i)]*pd_samples[i] + l_qb_val[string(i)]*qd_samples[i]  for i in load_buses)
                @show pval_approx

                perr_pos_mc[to_approx["quantity_index"]] = max(perr_pos_mc[to_approx["quantity_index"]],pval - pval_approx)
                perr_neg_mc[to_approx["quantity_index"]] = max(perr_neg_mc[to_approx["quantity_index"]],pval_approx - pval)

            elseif to_approx["quantity"] == "line_reactive_power"
                qvar = pm.var[:nw][0][:q][to_approx["quantity_index"]]
                @show qval = getvalue(qvar)

                linearization_coefficients = linearization_coefficients_list[approx_num]

                l0_val = linearization_coefficients["l0"]
                l_pb_val = linearization_coefficients["l_pb"]
                l_qb_val = linearization_coefficients["l_qb"]

                qval_approx = l0_val + sum(l_pb_val[string(i)]*(sum(pg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))  +  l_qb_val[string(i)]*(sum(qg_samples[j] for j in pm.ref[:nw][0][:bus_gens][i]))   for i in gen_buses) - sum(l_pb_val[string(i)]*pd_samples[i] + l_qb_val[string(i)]*qd_samples[i]  for i in load_buses)
                @show qval_approx

                qerr_pos_mc[to_approx["quantity_index"]] = max(qerr_pos_mc[to_approx["quantity_index"]],qval - qval_approx)
                qerr_neg_mc[to_approx["quantity_index"]] = max(qerr_neg_mc[to_approx["quantity_index"]],qval_approx - qval)
            end
        end


    end
    # end of monte_carlo

    approximation_errors = Dict{Int,Any}()

    for (ind,to_approx) in to_approx_list
        if to_approx["quantity"] == "line_real_power"
            errors = Dict{String,Float64}()
            errors["positive_error"] = perr_pos_mc[to_approx["quantity_index"]]
            errors["negative_error"] = perr_neg_mc[to_approx["quantity_index"]]
            errors["maximum_error"] = max(perr_pos_mc[to_approx["quantity_index"]],perr_neg_mc[to_approx["quantity_index"]])
            approximation_errors[ind] = errors
        elseif to_approx["quantity"] == "line_reactive_power"
            errors = Dict{String,Float64}()
            errors["positive_error"] = qerr_pos_mc[to_approx["quantity_index"]]
            errors["negative_error"] = qerr_neg_mc[to_approx["quantity_index"]]
            errors["maximum_error"] = max(qerr_pos_mc[to_approx["quantity_index"]],qerr_neg_mc[to_approx["quantity_index"]])
            approximation_errors[ind] = errors
        else println("Quantity not supported")
        end
    end

    return approximation_errors


end
