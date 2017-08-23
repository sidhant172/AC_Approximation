# SUPPORT FUNCTIONS
function get_current_solution(solution::Dict{String,Any}, model::JuMP.Model, to_approx, ind_gen, ind_bus, ind_branch)
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


    quantity = to_approx["quantity"]
    index = to_approx["quantity_index"]
    val = 0
    if quantity == "line_real_power"
        val_JuMP_var = getindex(model, :p)
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "line_reactive_power"
        val_JuMP_var = getindex(model, :q)
        val = getvalue(val_JuMP_var[index])
    elseif quantity == "bus_voltage_magnitude"
        val = vm_curr[index]
    else println("Approximating ", quantity, " is not supported.")
    end


    current_sol["pg"] = pg_curr
    current_sol["qg"] = qg_curr
    current_sol["pd"] = pd_curr
    current_sol["qd"] = qd_curr
    current_sol["p"] = p_curr
    current_sol["q"] = q_curr
    current_sol["vm"] = vm_curr
    current_sol["va"] = va_curr
    current_sol["val"] = val

    return current_sol
end     # end of get_current_solution




# function to create useful structure in network data
