using GLPKMathProgInterface
using PowerModels
using Distributions


max_iter = 100

solver = GLPKSolverLP()

network_data = PowerModels.parse_file("../pglib_opf_case240_pserc.m")

num_bus = length(network_data["bus"])
# result = run_dc_opf(network_data, solver)

result_dict = Dict{String,Any}()


dist = Normal(0,0.03)

perturbations = rand(dist,num_bus)

for iter=1:max_iter
    network = deepcopy(network_data)
    for (i,bus) in network["bus"]
        bus["pd"] = bus["pd"]*(1+rand(dist))
    end
    result_dict[string(iter)] = run_dc_opf(network,solver)
end


# # m = Model(solver=IpoptSolver(mu_init=1e-9))
# m = Model(solver=IpoptSolver())
# #
# @variable(m,x,start=sqrt(1/3))
# # @variable(m,x,start=0.9944097250233357)
# @variable(m,y,start=sqrt(2/3))
# # @variable(m,y,start=0.10559026622375049)
# @constraint(m,x^2+y^2<=1)
# @constraint(m,x+y>=1.1)
# @objective(m,Min,y+0.2*x)
# #
# status = solve(m)

# m = Model(solver=IpoptSolver(linear_solver="ma97"))




# network_data = PowerModels.parse_file("nesta_case300_ieee.m")
#
# # solver = IpoptSolver()
# solver = IpoptSolver(linear_solver="ma97")
#
# run_ac_opf(network_data,solver)
