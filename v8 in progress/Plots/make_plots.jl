using Plots
gr()

using MAT
using PowerModels


# filename = "nesta_case57_ieee.m"
# filename = "nesta_case30_as.m"
# filename = "nesta_case14_ieee.m"
filename = "case24_ieee_rts.m"
# filename = ARGS[1]


network_data = PowerModels.parse_file("../"string(filename))
num_branch = length(network_data["branch"])
num_bus = length(network_data["bus"])

inflation_factors = [0.1, 0.2, 0.3, 0.4]

# inflation_factors = [0.1]

vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_real0.05.mat")
# vars_jac = matread("../ptdf_matrices.mat")

pp = vars["coeff_p"]
active_buses = [i for i=1:num_bus if norm(pp[:,i])>1e-5]

# pp_jac = vars_jac["Hac_f"][1:38,1:24]
# pq_jac = vars_jac["Hac_f"][1:38,25:48]
# qp_jac = vars_jac["Hac_f"][39:76,1:24]
# qq_jac = vars_jac["Hac_f"][39:76,25:48]

# p = savefig(histogram2d(randn(10000),randn(10000),nbins=100),"mylot.pdf")

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
    # println("showing results for radius = ",inflation)
    # @show norm(pp[:,active_buses]-pp_jac[:,active_buses])
    # @show norm(pq[:,active_buses]-pq_jac[:,active_buses])
    # @show norm(qp[:,active_buses]-qp_jac[:,active_buses])
    # @show norm(qq[:,active_buses]-qq_jac[:,active_buses])
    # @show norm(p_err)
    # @show norm(q_err)
    # savefig(heatmap(active_buses,1:38,abs(pp[:,active_buses]-pp_jac[:,active_buses]),aspect_ratio=1),"heatmap_pp_"string(inflation)".pdf")
    # savefig(heatmap(active_buses,1:38,abs(pq[:,active_buses]-pq_jac[:,active_buses]),aspect_ratio=1),"heatmap_pq_"string(inflation)".pdf")
    # savefig(heatmap(active_buses,1:38,abs(qp[:,active_buses]-qp_jac[:,active_buses]),aspect_ratio=1),"heatmap_qp_"string(inflation)".pdf")
    # savefig(heatmap(active_buses,1:38,abs(qq[:,active_buses]-qq_jac[:,active_buses]),aspect_ratio=1),"heatmap_qq_"string(inflation)".pdf")
end


p_error_old = zeros(num_branch,length(inflation_factors))
q_error_old = zeros(num_branch,length(inflation_factors))
ctr = 0
olddirname = "/Users/SidhantMisra/Dropbox/Work/GitHub/AC_Approximation/v7"
for inflation in inflation_factors
    ctr = ctr+1
    vars = matread(string(olddirname)"/results"string(num_bus)"/matrix_forms/linear_approximations_real"string(inflation)".mat")
    pp = vars["coeff_p"]
    pq = vars["coeff_q"]
    p_err_old = vars["approx_error"]
    p_error_old[:,ctr] = p_err_old
    vars = matread(string(olddirname)"/results"string(num_bus)"/matrix_forms/linear_approximations_reactive"string(inflation)".mat")
    qp = vars["coeff_p"]
    qq = vars["coeff_q"]
    q_err_old = vars["approx_error"]
    q_error_old[:,ctr] = q_err_old
end





for i = 1:num_branch
    realfig = plot(inflation_factors,p_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",label="cg")
    plot!(inflation_factors,p_error_old[i,:],color=:red,label="gd")
    savefig(realfig,"plots"string(num_bus)"/approximation_error/approx_error_real_line_"string(i)".pdf")

    reactivefig = plot(inflation_factors,q_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",label="cg")
    plot!(inflation_factors,q_error_old[i,:],color=:red,label="gd")
    savefig(realfig,"plots"string(num_bus)"/approximation_error/approx_error_reactive_line_"string(i)".pdf")

    # savefig(plot(inflation_factors,p_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_real_line_"string(i)".pdf")
    # savefig(plot(inflation_factors,q_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_reactive_line_"string(i)".pdf")
end
