using JuMP
using Ipopt

function f(x...)
    return sum((x[i]-1)^2 for i in 1:length(x))
end

function fprime(g,x...)
    for i=1:length(x)
        g[i] = 2*(x[i]-1)
    end
end

# function fdprime(h,x...)
#     for i=1:length
# end

m = Model(solver=IpoptSolver(hessian_approximation="limited-memory"))

JuMP.register(m,:f,2,f,fprime)
# JuMP.register(m,:f,4,f,autodiff=true)

@variable(m,x[1:2])
# @variable(m,y)
# @variable(m,z)

# args = (x[1],x[2])
@NLobjective(m, Min, f(x[1],x[2]))

# @NLconstraint(m, f(x...) <= z)

status = solve(m)


# x = [3,2]
#
#
# @show f(x...)

# function f(x...)
#     return sum(x[i] for i in 1:length(x))
# end
#
# x = [1,2,3,4]
# arg = (x[1],x[2],x[3])
# # @show f(x[1],x[2],x[3])
# @show f(x...)
# function f(x; kwargs...)
#     @show c = kwargs[1][2]
#     return x + c
# end
#
#
# @show f(4; :c=>2)
