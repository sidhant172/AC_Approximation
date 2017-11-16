using JuMP
using Ipopt
using Plots
gr()


# thetavals =  0:0.1:10
# corr = zeros(length(thetavals))
#
#
# for i=1:length(thetavals)
#
#     theta = thetavals[i]
#
#     m = Model(solver=IpoptSolver())
#
#     @variable(m,0<=s<=1)
#
#     @NLobjective(m,Max,-s*log(s)-(1-s)*log(1-s) + 0.5*theta*(1-2s)^2)
#
#     status = solve(m)
#
#     @show sval = getvalue(s)
#
#     @show corr[i] = (exp(2*theta*(1-2*sval))-1)/(exp(2*theta*(1-2*sval))+1)
# end


c = 0.1


# for i=10:1000
i = 10000
    A = (1-c/i)*eye(i) + (c/i)*ones(i,i)
    @show minimum(eigvals(A))
# end



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
