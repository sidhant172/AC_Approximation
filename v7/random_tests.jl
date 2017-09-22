using JuMP
using Ipopt
using PowerModels

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
m = Model(solver=IpoptSolver())

@variable(m,-3<=x<=3,start=0.9)
@variable(m,-3<=y<=3,start=0.3)
# @variable(m,x,start=sqrt(1/3))
# @variable(m,y,start=sqrt(2/3))

@variable(m,l1>=0)
@variable(m,l2>=0)

@objective(m,Min,y)

# complementary slackness
@NLconstraint(m,l1*(x^2+y^2-1) == 0)
@NLconstraint(m,l2*(1.1-x-y) == 0)


# primal feasibility
@constraint(m,x^2+y^2<=1)
@constraint(m,x+y>=1.1)
# @NLconstraint(m,-1e-2 <= l1*(x^2+y^2-1) <=1e-2)
# @NLconstraint(m,-1e-2 <= l2*(1.1-x-y) <=1e-2)

# derivative conditions
@NLconstraint(m,0.2+l1*(2*x)+l2*(-1) == 0)

status=solve(m)




# network_data = PowerModels.parse_file("nesta_case300_ieee.m")
#
# # solver = IpoptSolver()
# solver = IpoptSolver(linear_solver="ma97")
#
# run_ac_opf(network_data,solver)
