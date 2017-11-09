using JuMP
using Ipopt

function f(x...)
    return (x[1]-1)^2 +x[2]^2
end

function fprime(g,x...)
    g[1]=2*(x[1]-1)
    g[2]=2*x[2]
end



m = Model(solver=IpoptSolver(hessian_approximation="limited-memory"))

JuMP.register(m,:f,2,f,fprime)

@variable(m,x[1:2])

args = (x[1],x[2])

@NLobjective(m, Min, f(x[1],x[2]) )

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
