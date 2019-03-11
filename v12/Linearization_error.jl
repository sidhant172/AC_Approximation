using PowerModels
using Ipopt
using JLD
using JuMP
using MAT


# eps = parse(Float64,ARGS[1])
# num_samples = parse(Int64,ARGS[2])
# # eps = 0.0
# # num_samples = 1
# # eps = 0.7
# # num_samples = 1000
#
# @show eps
# @show num_samples


solver = IpoptSolver(print_level=0)
# data = PowerModels.parse_file("Cases/case24.m")
# data = PowerModels.parse_file("Cases/case14.m")
data = PowerModels.parse_file("Cases/case9.m")
# data = PowerModels.parse_file("Cases/case3.m")


# data = PowerModels.parse_file("Cases/nesta_case3_lmbd.m")
# load_buses = [2,3]

load_buses = [load["load_bus"] for (i,load) in data["load"]]
gen_buses = [gen["gen_bus"] for (i,gen) in data["gen"]]
num_bus = length(data["bus"])
load_buses = setdiff(load_buses,gen_buses)

pm_opf = PowerModels.build_generic_model(data,ACPPowerModel,PowerModels.post_opf)
result_opf = PowerModels.solve_generic_model(pm_opf,solver)

data_opf = deepcopy(data)
for (i,gen) in data_opf["gen"]
    gen["pg"] = result_opf["solution"]["gen"][i]["pg"]
    gen["qg"] = result_opf["solution"]["gen"][i]["qg"]
end
for (i,bus) in data_opf["bus"]
    bus["vm"] = result_opf["solution"]["bus"][i]["vm"]
end

bus_p = Dict{Int64,Float64}()
bus_q = Dict{Int64,Float64}()
bus_pd = Dict{Int64,Float64}()
bus_qd = Dict{Int64,Float64}()
bus_pg = Dict{Int64,Float64}()
bus_qg = Dict{Int64,Float64}()

for l in load_buses
    bus_p[l], bus_q[l] = 0,0
    bus_pd[l], bus_qd[l] = 0,0
    bus_pg[l], bus_qg[l] = 0,0

    for (i,gen) in data["gen"]
        if gen["gen_bus"] == l
            bus_pg[l] += result_opf["solution"]["gen"][i]["pg"]
            bus_qg[l] += result_opf["solution"]["gen"][i]["qg"]
        end
    end
    for (i,load) in data["load"]
        if load["load_bus"] == l
            bus_pd[l] += load["pd"]
            bus_qd[l] += load["qd"]
        end
    end
    bus_p[l] = bus_pg[l] - bus_pd[l]
    bus_q[l] = bus_qg[l] - bus_qd[l]
end


bus_loads = Dict{Int64,Any}()
for l in load_buses
    bus_loads[l] = []
    for (i,load) in data["load"]
        if load["load_bus"] == l
            push!(bus_loads[l],i)
        end
    end
end

bus_gens = Dict{Int64,Any}()
for l in load_buses
    bus_gens[l] = []
    for (i,gen) in data["gen"]
        if gen["gen_bus"] == l
            push!(bus_gens[l],i)
        end
    end
end

inflation_factors = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]

num_branch = length(data["branch"])
num_bus = length(data["bus"])

p_error = zeros(num_branch)
q_error = zeros(num_branch)


for inflation_num=1:length(inflation_factors)
# for num = 1:num_samples
    @show num
    data_mod = deepcopy(data_opf)
    r1 = Dict{Int64,Any}()
    r2 = Dict{Int64,Any}()
    for l in load_buses
        r1[l] = 2*rand()-1
        r2[l] = 2*rand()-1
    end

    for l in load_buses
        for j in bus_loads[l]
            # @show eps*r1[l]*bus_pd[l]
            # @show eps*r2[l]*bus_qd[l]
            data_mod["load"][j]["pd"] = data_mod["load"][j]["pd"] + eps*r1[l]*bus_pd[l]
            data_mod["load"][j]["qd"] = data_mod["load"][j]["qd"] + eps*r2[l]*bus_qd[l]
            break
        end
        # data_mod["load"][bus_loads[l]]["pd"] += bus_pd[l]*r1[l]*eps
        # data_mod["load"][bus_loads[l]]["qd"] += bus_qd[l]*r2[l]*eps
    end

    pm = PowerModels.build_generic_model(data_mod, ACPPowerModel,PowerModels.post_pf)
    result = PowerModels.solve_generic_model(pm,solver)

    # filename = string("Cases/case24_ptdf_Vfixed.mat")
    # filename = string("Cases/case14_ptdf_Vfixed.mat")
    filename = string("Cases/case9_ptdf.mat")

    # filename = string("Cases/case3_lmbd_ptdf.mat")
    vars = matread(filename)
    jac = vars["Hac_f"]

    # for line = 1:num_branch
for line = 1:1
        f_bus,t_bus = data_mod["branch"][string(line)]["f_bus"],data_mod["branch"][string(line)]["t_bus"]

        # lp = jac[line,1:num_bus]
        # lq = jac[line,num_bus+1:end]
        # no_load_buses = setdiff(1:num_bus,load_buses)
        # lp[no_load_buses] = 0
        # lq[no_load_buses] = 0
        # p_actual_opf = getvalue(PowerModels.var(pm_opf,:p)[(line,f_bus,t_bus)])
        # l0 = p_actual_opf - sum(bus_p[l]*lp[l] + bus_q[l]*lq[l] for l in load_buses)

        p_actual = getvalue(PowerModels.var(pm,:p)[(line,f_bus,t_bus)])
        # lp = jac[line,1:num_bus]
        # lq = jac[line,num_bus+1:end]
        # no_load_buses = setdiff(1:num_bus,load_buses)
        # lp[no_load_buses] = 0
        # lq[no_load_buses] = 0
        # p_actual_opf = getvalue(PowerModels.var(pm_opf,:p)[(line,f_bus,t_bus)])
        # l0 = p_actual_opf - sum(bus_p[l]*lp[l] + bus_q[l]*lq[l] for l in load_buses)
        # p_predicted = l0 + sum(lp[l]*(bus_pg[l] - sum(data_mod["load"][j]["pd"] for j in bus_loads[l])) for l in load_buses) + sum(lq[l]*(bus_qg[l] - sum(data_mod["load"][j]["qd"] for j in bus_loads[l])) for l in load_buses)
        # p_error[line] += (p_predicted-p_actual)^2

        @show q_actual = getvalue(PowerModels.var(pm,:q)[(line,f_bus,t_bus)])
        # lp = jac[line+num_branch,1:num_bus]
        # lq = jac[line+num_branch,num_bus+1:end]
        # no_load_buses = setdiff(1:num_bus,load_buses)
        # lp[no_load_buses] = 0
        # lq[no_load_buses] = 0
        # q_actual_opf = getvalue(PowerModels.var(pm_opf,:q)[(line,f_bus,t_bus)])
        # l0 = q_actual_opf - sum(bus_p[l]*lp[l] + bus_q[l]*lq[l] for l in load_buses)
        # q_predicted = l0 + sum(lp[l]*(bus_pg[l] - sum(data_mod["load"][j]["pd"] for j in bus_loads[l])) for l in load_buses) + sum(lq[l]*(bus_qg[l] - sum(data_mod["load"][j]["qd"] for j in bus_loads[l])) for l in load_buses)
        # q_error[line] += (q_predicted-q_actual)^2






        # the commented code generates results for the optimal approximation obtained by PCE
        f_bus,t_bus = data_mod["branch"][string(line)]["f_bus"],data_mod["branch"][string(line)]["t_bus"]

        # active power approximation
        filename = string("Results_",string(num_bus),"/eps_",string(eps),"/active_line_",string(line),".jld")
        # filename = string("Results/eps_",string(0.01),"/active_line_",string(line),".jld")
        vars = JLD.load(filename)
        l0 = vars["l0"]
        lp = vars["lp"]
        lq = vars["lq"]
        p_predicted = l0 + sum(lp[l]*(bus_pg[l] - sum(data_mod["load"][j]["pd"] for j in bus_loads[l])) for l in load_buses) + sum(lq[l]*(bus_qg[l] - sum(data_mod["load"][j]["qd"] for j in bus_loads[l])) for l in load_buses)
        p_actual = getvalue(PowerModels.var(pm,:p)[(line,f_bus,t_bus)])
        p_error[line] += (p_predicted-p_actual)^2

        # reactive power approximation
        filename = string("Results_",string(num_bus),"/eps_",string(eps),"/reactive_line_",string(line),".jld")
        # filename = string("Results/eps_",string(0.01),"/reactive_line_",string(line),".jld")
        vars = JLD.load(filename)
        l0 = vars["l0"]
        lp = vars["lp"]
        lq = vars["lq"]
        q_predicted = l0 + sum(lp[l]*(bus_pg[l] - sum(data_mod["load"][j]["pd"] for j in bus_loads[l])) for l in load_buses) + sum(lq[l]*(bus_qg[l] - sum(data_mod["load"][j]["qd"] for j in bus_loads[l])) for l in load_buses)

        q_actual = getvalue(PowerModels.var(pm,:q)[(line,f_bus,t_bus)])
        q_error[line] += (q_predicted-q_actual)^2
        # end
    end

end # end of num_samples

p_error = p_error/num_samples
q_error = q_error/num_samples

# approximation error for optimal approximation from PCE
JLD.save(string("Results_",string(num_bus),"/eps_",string(eps),"/MC_error.jld"),"p_error",p_error,"q_error",q_error)

# approximation error for jacobian
# JLD.save(string("Results_",string(num_bus),"/eps_",string(eps),"/MC_error_jacobian.jld"),"p_error",p_error,"q_error",q_error)


println("JLD file written ...")
