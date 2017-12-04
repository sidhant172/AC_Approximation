using JuMP
using Ipopt
# using Plots
# gr()

n = 3

step = 0.05

minval = 100
thetamin = eye(n)
varmin = eye(n)

# for t12 = -gamma:0.02:gamma, t13 = -gamma:0.02:gamma, t14 = -gamma:0.02:gamma, t23 = -gamma:0.02:gamma, t24 = -gamma:0.02:gamma, t34 = -gamma:0.02:gamma
#         if (abs(t12)+abs(t13)+abs(t14) > gamma)
#             continue
#         elseif (abs(t12)+abs(t23)+abs(t24) > gamma)
#             continue
#         elseif (abs(t13)+abs(t23)+abs(t34) > gamma)
#             continue
#         elseif (abs(t14)+abs(t24)+abs(t34) > gamma)
#             continue
#         end
#
#         theta = eye(n)
#         theta[1,2]=t12
#         theta[1,3]=t13
#         theta[1,4]=t14
#         theta[2,3]=t23
#         theta[2,4]=t24
#         theta[3,4]=t34
#
#         for i=1:n,j=i+1:n
#             theta[j,i] = theta[i,j]
#         end
#
#
#         expval = zeros(2^n)
#         for l = 0:2^n-1
#             expval[l+1] = exp(sum(theta[i,j]*(2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1) for i=1:n for j=i+1:n))
#         end
#         z = sum(expval)
#
#         varmat = eye(n)
#         for i=1:n, j=i+1:n
#             varmat[i,j] = (1/z)*sum( (2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1)*expval[l+1] for l=0:2^n-1 )
#             varmat[j,i] = varmat[i,j]
#         end
#
#         smallest_eig = minimum(eigvals(varmat))
#
#         if minval > smallest_eig
#             minval = smallest_eig
#             thetamin = theta
#             varmin = varmat
#         end
#
# end


d = Dict{Array{Float64,1},Float64}()

for t23 = -gamma:step:gamma
    for t12 = -(gamma-abs(t23)):step:(gamma-abs(t23)), t13 = -(gamma-abs(t23)):step:(gamma-abs(t23))
        if (abs(t13)+abs(t12) > gamma)
            continue
        end
        theta = eye(n)
        theta[2,3]=t23
        theta[3,2] = theta[2,3]
        theta[1,3]=t13
        theta[3,1] = theta[1,3]
        theta[1,2]=t12
        theta[2,1] = theta[1,2]


        expval = zeros(2^n)
        for l = 0:2^n-1
            expval[l+1] = exp(sum(theta[i,j]*(2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1) for i=1:n for j=i+1:n))
        end
        z = sum(expval)

        varmat = eye(n)
        for i=1:n, j=i+1:n
            varmat[i,j] = (1/z)*sum( (2*digits(l,2,n)[i]-1)*(2*digits(l,2,n)[j]-1)*expval[l+1] for l=0:2^n-1 )
            varmat[j,i] = varmat[i,j]
        end

        smallest_eig = minimum(eigvals(varmat))

        d[[t23,t12,t13]] = smallest_eig

        if minval > smallest_eig
            minval = smallest_eig
            thetamin = theta
            varmin = varmat
        end

    end
end




println("Minimum eigenvalue = ",minval)
println(thetamin)
println(varmin)
println("1-tanh(gamma)=",1 - tanh(gamma))
