using JuMP
using Ipopt

function create_model()
    m = Model(solver = IpoptSolver())
    @variable(m,x>=0)
    @objective(m,Min,x)
    status = solve(m)
    return m
end

function get_sol(m)
    println("Printing x value = ", getvalue(x))
end
