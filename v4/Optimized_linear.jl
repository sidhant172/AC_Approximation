using PowerModels
using JuMP
using Ipopt
using Clp

include("opf_mod.jl")
include("support_functions.jl")



################################################################################
# OPTIONS

# algorithm parameters
cnst_gen_max_iter  = 1000   # max iterations for constraint generation
tol = 1e-3   # convergence tolerance

# operational conditions
gen_inflation = 0.1     # defining range of loading conditions
load_inflation = 0.1    # defining range of generation conditions


# choose transmission lines and buses to find best linearazion over
lines = [(18,11,13)]

buses = [2]

quantity_to_approx = "line_real_power"
# quantity_to_approx = "line_reactive_power"
# quantity_to_approx = "bus_voltage_magnitude"

################################################################################





network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("case118.m")


line = (18,11,13)
# line = (11,5,11)


solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
# solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")

network_data["line"] = line

# Extracting useful structures
branch = network_data["branch"]
bus = network_data["bus"]
gen = network_data["gen"]


######### Defining indices and parameters ######################################
ind_bus = [parse(Int,key) for (key,b) in bus]
ind_branch = [parse(Int,key) for (key,b) in branch]
ind_gen = [parse(Int,key) for (key,b) in gen]

n_bus = length(ind_bus)
n_branch = length(ind_branch)
n_gen = length(ind_gen)

################################################################################


############## Mapping generators to buses #####################################
gen_buses = [gen[i]["gen_bus"] for i in keys(gen)]
gen_buses = unique(gen_buses)

gens_at_bus = Dict{Int64,Array{Int,1}}()

load_buses = [parse(Int64,i) for i in keys(bus) if abs(bus[i]["pd"]) + abs(bus[i]["qd"]) > 1e-2]

active_buses = union(gen_buses,load_buses)

for i in gen_buses
    gens_at_bus[i] = [gen[j]["index"] for j in keys(gen) if gen[j]["gen_bus"] == i]
end

################################################################################



############# vary generators around nominal opf solution ######################
output = run_ac_opf(network_data, solver_ipopt)

pg_init = Dict{Int,Float64}()
qg_init = Dict{Int,Float64}()
baseMVA = output["solution"]["baseMVA"]

for i in ind_gen
    pg_init[i] = output["solution"]["gen"][string(i)]["pg"]/baseMVA
    qg_init[i] = output["solution"]["gen"][string(i)]["qg"]/baseMVA
    network_data["gen"][string(i)]["pmax"] = max(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
    network_data["gen"][string(i)]["pmin"] = min(pg_init[i]*(1+gen_inflation),pg_init[i]*(1-gen_inflation))
    network_data["gen"][string(i)]["qmax"] = max(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
    network_data["gen"][string(i)]["qmin"] = min(qg_init[i]*(1+gen_inflation),qg_init[i]*(1-gen_inflation))
end
################################################################################



# tic()   # start the clock

m = Model(solver=ClpSolver())   # JuMP model for master problem

@variable(m, z >=0, start=0)  # maximum error

# linearization coefficients (currently for one line)
@variable(m,l_pb[i in active_buses], start = 0)
@variable(m,l_qb[i in active_buses], start = 0)
@variable(m, l0, start = 0) # constant term


@objective(m,Min,z)


l0_val = 0
l_pb_val = Dict{String,Float64}()
l_qb_val = Dict{String,Float64}()
for i in active_buses
    l_pb_val[string(i)] = 0.0
    l_qb_val[string(i)] = 0.0
end

converged_flag = 0
pos = 0


lin_coeff = Dict{String, Any}()     # container for storing current value of linearization coefficients

#lp_err = Float64[];
#nlp_err = Float64[];

# details about which quantity and are running the constraint generation for
to_approx = Dict{String,Any}()
to_approx["quantity"] = quantity_to_approx
to_approx["index"] = line
# network_data["to_approx"] = to_approx   # add to network_data
network_data["quantity"] = quantity_to_approx
network_data["quantity_index"] = line



for iter=1:cnst_gen_max_iter   # total iterations of constraint generation
    status = solve(m)
    z_val  = getvalue(z)
    l0_val = getvalue(l0)
    for i in active_buses
        l_pb_val[string(i)] = getvalue(l_pb[i])
        l_qb_val[string(i)] = getvalue(l_qb[i])
    end

    # add linearization coefficients to network_data
    # lin_coeff["l0"] = l0_val
    # lin_coeff["l_pb"] = l_pb_val
    # lin_coeff["l_qb"] = l_qb_val
    # network_data["lin_coeff"] = lin_coeff
    network_data["l0"] = l0_val
    network_data["l_pb"] = l_pb_val
    network_data["l_qb"] = l_qb_val

    # l_pb_val = getvalue(l_pb)
    # l_qb_val = getvalue(l_qb)
    #push!(lp_err,abs(z_val))
    #push!(nlp_err,0)

    converged_flag = 1

################################################################################
    network_data["direction"] = 0   # direction of maximization
    (result, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)

    @show z_val - result["objective"]
    #nlp_err[iter] = abs(z_val - s["objective"])

    if (z_val - result["objective"]) < -tol
        converged_flag = 0

        current_sol = get_current_solution(result["solution"], pm_model, to_approx)

        @constraint(m, current_sol["val"] - (l0 +
            sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[i]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[i]) ) for i in gen_buses)
        -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
            <= z)
    end
################################################################################

################################################################################
    network_data["direction"] = 1   # direction of maximization
    (result, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)

    @show z_val - result["objective"]
    #nlp_err[iter] = abs(z_val - s["objective"])

    if (z_val - result["objective"]) < -tol
        converged_flag = 0

        current_sol = get_current_solution(result["solution"], pm_model, to_approx)

        @constraint(m, -current_sol["val"] + (l0 +
        sum(l_pb[i]*( sum(current_sol["pg"][j] for j in gens_at_bus[i]) ) + l_qb[i]*( sum(current_sol["qg"][j] for j in gens_at_bus[i]) ) for i in gen_buses)
    -   sum(l_pb[i]*current_sol["pd"][i]  + l_qb[i]*current_sol["qd"][i] for i in load_buses) )
        <= z)
    end
################################################################################

@show iter
@show z_val

    if converged_flag == 1
        println("Constraint generation converged!")
        break
    end
end     # end of constraint generation loop
################################################################################
# time = toc()
# @show time  # stop the clock
