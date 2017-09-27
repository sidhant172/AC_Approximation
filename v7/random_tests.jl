

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
println(ARGS[1])
println(ARGS[2])
@show ARGS[1]


@show parse(Int64,ARGS[1]) + parse(Int64,ARGS[2])



# network_data = PowerModels.parse_file("nesta_case300_ieee.m")
#
# # solver = IpoptSolver()
# solver = IpoptSolver(linear_solver="ma97")
#
# run_ac_opf(network_data,solver)
