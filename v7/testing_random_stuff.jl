using JuMP
using Ipopt
using GLPKMathProgInterface

eps = 1e-3
tol = 1e-8

m = Model(solver = GLPKSolverLP())
@variable(m,z>=0)
@variable(m,-10<=a<=10)
@variable(m,-10<=b<=10)

# @constraint(m,a==1)

@objective(m,Min,z)


converged = 0

for ctr = 1:1000

    status = solve(m)

    aval = getvalue(a)
    bval = getvalue(b)
    z_val = getvalue(z)


    converged = 1

    ml = Model(solver=IpoptSolver(print_level=0))
    @variable(ml,1-eps<=xl<=1+eps)
    @NLobjective(ml,Max,xl^3-aval-bval*xl)
    status = solve(ml)

    xl_val = getvalue(xl)
    @show  z_val - (xl_val^3-aval-bval*xl_val)
    if z_val - (xl_val^3-aval-bval*xl_val) < -tol
        converged = 0
        @constraint(m,z >= xl_val^3 - a - b*xl_val)
    end

    mu = Model(solver=IpoptSolver(print_level=0))
    @variable(mu,1-eps<=xu<=1+eps)
    @NLobjective(mu,Max,-xu^3+aval+bval*xu)
    status = solve(mu)

    xu_val = getvalue(xu)
    @show z_val - (-xu_val^3+aval+bval*xu_val)
    if z_val - (-xu_val^3+aval+bval*xu_val) < -tol
        converged = 0
        @constraint(m,z >= a + b*xu_val - xu_val^3 )
    end

    if converged == 1
        println("Solution converged!")
        break
    end
end



mjl = Model(solver=IpoptSolver(print_level=0))
@variable(mjl,1-eps<=xjl<=1+eps)
@NLobjective(mjl,Max,xjl^3+2-3*xjl)
status = solve(mjl)
xjl_val = getvalue(xjl)
@show xjl_val^3+1-3*xjl_val

mju = Model(solver=IpoptSolver(print_level=0))
@variable(mju,1-eps<=xju<=1+eps)
@NLobjective(mju,Max,-xju^3-2+3*xju)
status = solve(mju)
xju_val = getvalue(xju)
@show -xju_val^3-2+3*xju_val
