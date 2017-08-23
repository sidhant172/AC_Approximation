using PowerModels
using JuMP
using Ipopt
using Clp
using MAT

include("opf_mod.jl")
include("support_functions.jl")


################################################################################
# Algorithm parameters
cnst_gen_max_iter  = 1000   # max iterations for constraint generation
tol = 1e-3   # convergence tolerance
gen_inflation = 0.2
load_inflation = 0.2
################################################################################



network_data = PowerModels.parse_file("case24_ieee_rts.m")
# network_data = PowerModels.parse_file("case118.m")


line = (18,11,13)
# line = (11,5,11)


solver_ipopt = IpoptSolver(print_level=0) # , linear_solver="ma97"
# solver_ipopt = IpoptSolver(print_level=0, linear_solver="ma97")



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



m = Model(solver=ClpSolver())   # JuMP model for master problem

@variable(m,z >=0, start=0)  # maximum error

# linearization coefficients (currently for one line)
@variable(m,l_pg[i in ind_gen])
@variable(m,l_qg[i in ind_gen])
@variable(m,l_pd[i in ind_bus])
@variable(m,l_qd[i in ind_bus])
@variable(m, l0) # constant term
# @variable(m,l_qd)

@objective(m,Min,z)

@constraint(m, l_pg[1] == 0)
@constraint(m, l_qg[1] == 0)

l_pd_val  = zeros(n_bus)
l_qd_val  = zeros(n_bus)
l_pg_val  = zeros(n_gen)
l_qg_val  = zeros(n_gen)
l0_val = 0

converged_flag = 0
pos = 0

lp_err = Float64[];
nlp_err = Float64[];

for iter=1:cnst_gen_max_iter   # total iterations of constraint generation

    status = solve(m)
    z_val  = getvalue(z)
    l0_val = getvalue(l0)
    l_pd_val = getvalue(l_pd)
    l_qd_val = getvalue(l_qd)
    l_pg_val = getvalue(l_pg)
    l_qg_val = getvalue(l_qg)

    push!(lp_err,abs(z_val))
    push!(nlp_err,0)

    converged_flag = 1


################################################################################
    pos = 0     # direction of maximization

    (s, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)

    solution = s["solution"]

    @show z_val - s["objective"]
    nlp_err[iter] = abs(z_val - s["objective"])

    if (z_val - s["objective"]) < -tol
        converged_flag = 0

        current_sol = get_current_solution(solution, pm_model)

        @constraint(m, current_sol["pval"] - (l0 + sum(l_pg[i]*current_sol["pg"][i] + l_qg[i]*current_sol["qg"][i] for i in ind_gen)
        + sum(l_pd[i]*current_sol["pd"][i] + l_qd[i]*current_sol["qd"][i] for i in ind_bus))  <= z)
    end

################################################################################

################################################################################
    pos = 1     # direction of maximization

    (s, pm_model) = run_ac_opf_mod(network_data,solver_ipopt)
    # s = run_ac_opf_mod(network_data,solver_ipopt)
    solution = s["solution"]
    # @show s

    nlp_err[iter] = max(nlp_err[iter],abs(z_val - s["objective"]))

    @show z_val - s["objective"]

    if (z_val - s["objective"]) < -tol
        converged_flag = 0

        current_sol = get_current_solution(solution, pm_model)

        @constraint(m, -current_sol["pval"] + (l0 + sum(l_pg[i]*current_sol["pg"][i] + l_qg[i]*current_sol["qg"][i] for i in ind_gen)
        + sum(l_pd[i]*current_sol["pd"][i] + l_qd[i]*current_sol["qd"][i] for i in ind_bus))  <= z)
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

# Output errors to a Matlab file
matwrite("err.mat", Dict(
	"lp_err" => lp_err,
	"nlp_err" => nlp_err
))
