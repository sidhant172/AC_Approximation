using JuMP
using Ipopt
# using Plots
# gr()

n = 6
gamma = 0.1
eps = 1
###### smallest eigen value
m = Model(solver=IpoptSolver(linear_solver="ma57",tol=1e-14))

@variable(m,-2<=theta[i=1:n,j=i+1:n]<=2,start = (1+eps*(2*rand()-1))*gamma/(n-1))
@variable(m,slack[i=1:n,j=i+1:n]>=0,start=(1+eps*(2*rand()-1))*gamma/(n-1))
@variable(m,x[1:n],start=1/sqrt(n))
# setvalue(x[1],1/sqrt(3))
# setvalue(x[2],-1/sqrt(3))
@variable(m,z>=0)
@variable(m,varmat[i=1:n,j=i+1:n])
@variable(m,expval[0:2^n-1]>=0)

@NLconstraint(m,expvaldef[l=0:2^n-1],expval[l] == exp(sum(theta[i,j]*(2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1) for i=1:n for j=i+1:n)))
@NLconstraint(m,z == sum(expval[l] for l=0:2^n-1))
@NLconstraint(m,varmatdef[i=1:n,j=i+1:n],varmat[i,j] == (1/z)*sum( (2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1)*expval[l] for l=0:2^n-1 ) )

@NLconstraint(m,sum(x[i]^2 for i=1:n) == 1)
@constraint(m,slackdefplus[i=1:n,j=i+1:n],theta[i,j] <= slack[i,j])
@constraint(m,slackdefminus[i=1:n,j=i+1:n],theta[i,j] >= -slack[i,j])
@constraint(m,l1norm[i=1:n],sum(slack[i,j] for j=i+1:n) + sum(slack[j,i] for j=1:i-1) <= gamma)

@NLobjective(m,Min,sum(x[i]*x[j]*varmat[i,j] for i=1:n for j=i+1:n))

status=solve(m)

@show 1  + 2*getobjectivevalue(m)
@show 1 - (tanh(gamma))
@show getvalue(x)
@show getvalue(theta)


mat = eye(n)
for i=1:n, j=i+1:n
    mat[i,j] = getvalue(varmat[i,j])
    mat[j,i] = mat[i,j]
end
@show mat



# function var(v,n)
#     theta =
#     zval = 0
#     expval = zeros(2^n)
#     for m = 0:2^n-1
#         b = digits(m,2,n)
#         expval[m] = sum(theta[i,j]*(2*b[i]-1)*(2*b[j]-1) for i=1:n for j=i+1:n)
#     end
#
#     Z = sum(expval) # partition function
#
#
#
#
#     sigmamat = zeros(n,n)
#     for i=1:n, j=i+1:n
#         sigmaij = 0
#         for m=0:2^n-1
#             b = digits(m,2,n)
#             sigmaij = sigmaij + (2*b[i]-1)*(2*b[j]-1)*expval[m]
#         end
#         sigmamat[i,j] = sigmaij/Z
#     end
#
#
#     return sum(x[i]^2 for i=1:n) + sum(sigmamat[i,j]*x[i]*x[j] for i=1:n for j=i+1:n)
# end
#
# thetavals =  0:0.01:5
# # corr = zeros(length(thetavals))
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
#     # @show corr[i] = (exp(2*theta*(1-2*sval))-1)/(exp(2*theta*(1-2*sval))+1)
#     @show corr[i] = tanh(theta*(1-2*sval))
# end


#
# # myplot = plot(thetavals,corr)
# myplot = plot(thetavals,(corr-1).*exp(2*thetavals))
# # *exp(0.5*thetavals)
# savefig(myplot,"corr_vs_theta.pdf")

# c = 0.1


# for i=10:1000
# # i = 10000
#     A = (1-c/(i-1))*eye(i) + (c/(i-1))*ones(i,i)
#     @show minimum(eigvals(A))
# # end



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
