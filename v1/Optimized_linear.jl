using PowerModels
using JuMP
using Ipopt
using Clp

cnst_gen_max_iter  = 1000

line = (18,11,13)
# line = (11,5,11)

# solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")

include("opf_mod.jl")

network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("case118.m")



# Extracting useful structures
branch = network_data["branch"]
bus = network_data["bus"]
gen = network_data["gen"]


######### Defining indices to avoid indexing problems due to removing elements in the data strcuture while creating contingencies #########
ind_bus = [parse(Int,key) for (key,b) in bus]
ind_branch = [parse(Int,key) for (key,b) in branch]
ind_gen = [parse(Int,key) for (key,b) in gen]

n_bus = length(ind_bus)
n_branch = length(ind_branch)
n_gen = length(ind_gen)

#######################################################
gen_buses = [dict["gen_bus"] for (i,dict) in gen]   # Find out which are generator buses, there are multiple per bus
gen_buses = unique(gen_buses)   # Remove duplicates due to multiple generators per node

bus_gen = Dict{Int,Array{Int,1}}()
for i in gen_buses
    buses = []
    for (gen_num,data) in gen
        if data["gen_bus"] == i
            buses = [buses;data["index"]]
        end
    end
    bus_gen[i] = buses
end
#######################################################

bus_arcs = Dict{Int,Array{Int,1}}()
for i in ind_bus
    arcs = []
    for (branch_num,data) in branch
        if data["f_bus"] == i
            arcs = [arcs;parse(Int64,branch_num)]
        end
        if data["t_bus"] == i
            arcs = [arcs;parse(Int64,branch_num)]
        end
    end
    bus_arcs[i] = arcs
end
#######################################################
bus_neighbors = Dict{Int,Array{Int,1}}()
for i in ind_bus
    neighbors = []
    for (branch_num,data) in branch
        if data["f_bus"] == i
            neighbors = [neighbors;data["t_bus"]]
        end
        if data["t_bus"] == i
            neighbors = [neighbors;data["f_bus"]]
        end
    end
    bus_neighbors[i] = neighbors
end
#######################################################




############# vary generators around nominal opf solution ############
output = run_ac_opf(network_data, solver_ipopt)

pg_init = Dict{Int,Float64}()
qg_init = Dict{Int,Float64}()
baseMVA = output["solution"]["baseMVA"]

for i in ind_gen
    pg_init[i] = output["solution"]["gen"][string(i)]["pg"]/baseMVA
    qg_init[i] = output["solution"]["gen"][string(i)]["qg"]/baseMVA
    network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*1.2,pg_init[i]*0.8)
    network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*1.2,pg_init[i]*0.8)
    network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*1.2,qg_init[i]*0.8)
    network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*1.2,qg_init[i]*0.8)
end
#####################################################################


# Find tight line constraints in the opf
# output = PowerModels.run_ac_opf(network_data, solver_ipopt)
# pm = build_generic_model(network_data, ACPPowerModel, PowerModels.post_opf)
# output = solve_generic_model(pm, solver_ipopt, solution_builder = PowerModels.get_solution)
#
# p = getindex(pm.model, :p)
# pval = getvalue(p[:])
#
# index_tight_line = []
# for i in ind_branch
#     if abs(pval[i]) - branch[string(i)][""]
# end






m = Model(solver=ClpSolver())

@variable(m,z >=0, start=0)  # maximum error

# linearization coefficients (currently for one line)
@variable(m,l_pg[i in ind_gen])
@variable(m,l_qg[i in ind_gen])
@variable(m,l_pd[i in ind_bus])
@variable(m,l_qd[i in ind_bus])
# @variable(m,l_qd)

@objective(m,Min,z)

l_pd_val  = zeros(n_bus)
l_qd_val  = zeros(n_bus)
l_pg_val  = zeros(n_gen)
l_qg_val  = zeros(n_gen)

converged_flag = 0
pos = 0

for iter=1:cnst_gen_max_iter   # total iterations of constraint generation

    status = solve(m)
    z_val  = getvalue(z)
    l_pd_val = getvalue(l_pd)
    l_qd_val = getvalue(l_qd)
    l_pg_val = getvalue(l_pg)
    l_qg_val = getvalue(l_qg)



    converged_flag = 1


################################################################################
    pos = 0     # direction of maximization

    (s, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)

    # s = run_ac_opf_mod(network_data,solver_ipopt)

    # @show getindex(pm_model, :p)

    solution = s["solution"]

    @show z_val - s["objective"]

    if (z_val - s["objective"]) < -1e-2
        converged_flag = 0

        baseMVA = s["solution"]["baseMVA"]

        pg_curr = Dict{Int,Float64}()
        qg_curr = Dict{Int,Float64}()
        vm_curr = Dict{Int,Float64}()
        va_curr = Dict{Int,Float64}()
        p_curr = Dict{Array{Int,1},Float64}()
        q_curr = Dict{Array{Int,1},Float64}()
        pd_curr = Dict{Int,Float64}()
        qd_curr = Dict{Int,Float64}()

        p = getindex(pm_model, :p)
        pval = getvalue(p[line])
        pd = getindex(pm_model, :pd)
        qd = getindex(pm_model, :pd)

        for i in ind_gen
            pg_curr[i] = solution["gen"][string(i)]["pg"]/baseMVA
            qg_curr[i] = solution["gen"][string(i)]["qg"]/baseMVA
        end

        for i in ind_bus
            vm_curr[i] = solution["bus"][string(i)]["vm"]
            pd_curr[i] = getvalue(pd[i])
            qd_curr[i] = getvalue(qd[i])
        end


        @constraint(m, pval - (sum(l_pg[i]*pg_curr[i] + l_qg[i]*qg_curr[i] for i in ind_gen)
        + sum(l_pd[i]*pd_curr[i] + l_qd[i]*qd_curr[i] for i in ind_bus))  <= z)

    end

################################################################################

################################################################################
    pos = 1     # direction of maximization

    (s, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)
    # s = run_ac_opf_mod(network_data,solver_ipopt)
    solution = s["solution"]
    # @show s

    @show z_val - s["objective"]

    if (z_val - s["objective"]) < -1e-2
        converged_flag = 0

        baseMVA = s["solution"]["baseMVA"]

        pg_curr = Dict{Int,Float64}()
        qg_curr = Dict{Int,Float64}()
        vm_curr = Dict{Int,Float64}()
        va_curr = Dict{Int,Float64}()
        p_curr = Dict{Array{Int,1},Float64}()
        q_curr = Dict{Array{Int,1},Float64}()
        pd_curr = Dict{Int,Float64}()
        qd_curr = Dict{Int,Float64}()

        p = getindex(pm_model, :p)
        pval = getvalue(p[line])
        pd = getindex(pm_model, :pd)
        qd = getindex(pm_model, :pd)

        for i in ind_gen
            pg_curr[i] = solution["gen"][string(i)]["pg"]/baseMVA
            qg_curr[i] = solution["gen"][string(i)]["qg"]/baseMVA
        end

        for i in ind_bus
            vm_curr[i] = solution["bus"][string(i)]["vm"]
            pd_curr[i] = getvalue(pd[i])
            qd_curr[i] = getvalue(qd[i])
        end


        @constraint(m, - pval + (sum(l_pg[i]*pg_curr[i] + l_qg[i]*qg_curr[i] for i in ind_gen)
        + sum(l_pd[i]*pd_curr[i] + l_qd[i]*qd_curr[i] for i in ind_bus))  <= z)

    end

################################################################################

@show iter
@show z_val

    if converged_flag == 1
        println("Constraint generation converged!")
        break
    end
end
