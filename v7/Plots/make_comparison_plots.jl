using Plots
gr()

using MAT
using PowerModels


# filename = "nesta_case57_ieee.m"
# filename = "nesta_case30_as.m"
# filename = "nesta_case14_ieee.m"
filename = "case24_ieee_rts.m"

network_data = PowerModels.parse_file("../"string(filename))
num_branch = length(network_data["branch"])
num_bus = length(network_data["bus"])

inflation_factors = [0.05, 0.1, 0.2, 0.3, 0.4]

# find active buses
vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_real0.05.mat")
pp = vars["coeff_p"]
active_buses = [i for i=1:num_bus if norm(pp[:,i])>1e-5]


################################################################################
# reading error of Taylor approximation
p_error_jac = zeros(num_branch,length(inflation_factors))
q_error_jac = zeros(num_branch,length(inflation_factors))

ctr = 0
for inflation in inflation_factors
    ctr = ctr+1
    vars = matread("../resultsJacobian/jacobian_error_real_"string(inflation)".mat")
    p_error_jac[:,ctr] = vars["err"]
    vars = matread("../resultsJacobian/jacobian_error_reactive_"string(inflation)".mat")
    q_error_jac[:,ctr] = vars["err"]
end
################################################################################


################################################################################
# reading error of optimal linear approximation
p_error = zeros(num_branch,length(inflation_factors))
q_error = zeros(num_branch,length(inflation_factors))

pps = Dict{String,Any}()
pqs = Dict{String,Any}()
qps = Dict{String,Any}()
qqs = Dict{String,Any}()


ctr = 0
for inflation in inflation_factors
    ctr = ctr+1
    vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_real"string(inflation)".mat")
    pp = vars["coeff_p"]
    pq = vars["coeff_q"]
    p_err = vars["approx_error"]
    p_error[:,ctr] = p_err
    vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_reactive"string(inflation)".mat")
    qp = vars["coeff_p"]
    qq = vars["coeff_q"]
    q_err = vars["approx_error"]
    q_error[:,ctr] = q_err

    # saving everyting
    pps[string(convert(Int64,100*inflation))] = pp
    pqs[string(convert(Int64,100*inflation))] = pq
    qps[string(convert(Int64,100*inflation))] = qp
    qqs[string(convert(Int64,100*inflation))] = qq
end
################################################################################



# for i = 1:num_branch
#     perr_plot = plot(inflation_factors,p_error_jac[i,:],linewidth = 2, xlabel="Radius",ylabel="Maximum approximation error (pu)",lab = "Taylor")
#     plot!(inflation_factors, p_error[i,:], linewidth = 2, lab = "Optimal")
#     savefig(perr_plot, "plotsJacobian/approximation_error_comparison/approx_error_real_line_"string(i)".pdf")
#
#     qerr_plot = plot(inflation_factors,q_error_jac[i,:], linewidth = 2, xlabel="Radius",ylabel="Maximum approximation error (pu)",lab = "Taylor")
#     plot!(inflation_factors, q_error[i,:], linewidth = 2, lab = "Optimal")
#     savefig(qerr_plot, "plotsJacobian/approximation_error_comparison/approx_error_reactive_line_"string(i)".pdf")
#     # savefig(plot(inflation_factors,p_error_jac[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plotsJacobian/approximation_error_comparison/approx_error_real_line_"string(i)".pdf")
#     # savefig(plot(inflation_factors,q_error_jac[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plotsJacobian/approximation_error_comparison/approx_error_reactive_line_"string(i)".pdf")
# end

# combined_plot = plot(size = (1400,900), layout=4, right_margin=20px, left_margin=20px, top_margin= 20px, bottom_margin=10px, [inflation_factors, inflation_factors, inflation_factors, inflation_factors],
#     [p_error_jac[1,:], q_error_jac[1,:], p_error_jac[18,:], q_error_jac[18,:]], linewidth=2, xlabel = "Radius", ylabel = "Maximum approximation error (pu)", lab="Taylor"
#     , annotations=[(0.1,0.03,text("Line #1 active",:left)) (0.1,0.24,text("Line #1 reactive",:left)) (0.1,0.19,text("Line #18 active",:left)) (0.1,0.15,text("Line #18 reactive",:left))] )
#
#     plot!(size = (1400,900), layout=4, right_margin=20px, left_margin=20px, top_margin= 20px,  bottom_margin=10px, [inflation_factors, inflation_factors, inflation_factors, inflation_factors],
#         [p_error[1,:], q_error[1,:], p_error[18,:], q_error[18,:]], linewidth=2, xlabel = "Radius", ylabel = "Maximum approximation error (pu)", lab="Optimal")
#
#     savefig(combined_plot, "plotsJacobian/approximation_error_comparison/combinedfigure.pdf")
