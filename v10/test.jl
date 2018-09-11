using JuMP
using KNITRO

solver = KnitroSolver()

m  = Model(solver=solver)
@variable(m,x>=0)

for i =1:10
    @objective(m,Min,x)
    status = solve(m)
end
