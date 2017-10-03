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
end
################################################################################

for i = 1:num_branch
    perr_plot = plot(inflation_factors,p_error_jac[i,:],linewidth = 2, xlabel="Radius",ylabel="Maximum approximation error (pu)",lab = "Taylor")
    plot!(inflation_factors, p_error[i,:], linewidth = 2, lab = "Optimal")
    savefig(perr_plot, "plotsJacobian/approximation_error_comparison/approx_error_real_line_"string(i)".pdf")

    qerr_plot = plot(inflation_factors,q_error_jac[i,:], linewidth = 2, xlabel="Radius",ylabel="Maximum approximation error (pu)",lab = "Taylor")
    plot!(inflation_factors, q_error[i,:], linewidth = 2, lab = "Optimal")
    savefig(qerr_plot, "plotsJacobian/approximation_error_comparison/approx_error_reactive_line_"string(i)".pdf")
    # savefig(plot(inflation_factors,p_error_jac[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plotsJacobian/approximation_error_comparison/approx_error_real_line_"string(i)".pdf")
    # savefig(plot(inflation_factors,q_error_jac[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plotsJacobian/approximation_error_comparison/approx_error_reactive_line_"string(i)".pdf")
end
