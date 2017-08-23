# SUPPORT FUNCTIONS
function get_current_solution(solution::Dict{String,Any}, model::JuMP.Model)
    baseMVA = solution["baseMVA"]

    current_sol = Dict{String,Any}()

    pg_curr = Dict{Int,Float64}()
    qg_curr = Dict{Int,Float64}()
    vm_curr = Dict{Int,Float64}()
    va_curr = Dict{Int,Float64}()
    p_curr = Dict{Array{Int,1},Float64}()
    q_curr = Dict{Array{Int,1},Float64}()
    pd_curr = Dict{Int,Float64}()
    qd_curr = Dict{Int,Float64}()

    p = getindex(model, :p)
    pval = getvalue(p[line])

    pd = getindex(model, :pd)
    qd = getindex(model, :qd)

    for i in ind_gen
        pg_curr[i] = solution["gen"][string(i)]["pg"]/baseMVA
        qg_curr[i] = solution["gen"][string(i)]["qg"]/baseMVA
    end

    for i in ind_bus
        vm_curr[i] = solution["bus"][string(i)]["vm"]
        pd_curr[i] = getvalue(pd[i])
        qd_curr[i] = getvalue(qd[i])
    end

    current_sol["pg"] = pg_curr
    current_sol["qg"] = qg_curr
    current_sol["pd"] = pd_curr
    current_sol["qd"] = qd_curr
    current_sol["p"] = p_curr
    current_sol["q"] = q_curr
    current_sol["vm"] = vm_curr
    current_sol["va"] = va_curr
    current_sol["pval"] = pval

    return current_sol
end     # end of get_current_solution
