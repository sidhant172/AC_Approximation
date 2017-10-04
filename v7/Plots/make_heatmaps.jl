using Plots
gr()

using MAT
using PowerModels


# filename = "nesta_case57_ieee.m"
# filename = "nesta_case30_as.m"
# filename = "nesta_case14_ieee.m"

# filename = ARGS[1]

filename = "case24_ieee_rts.m"

network_data = PowerModels.parse_file("../"string(filename))
num_branch = length(network_data["branch"])
num_bus = length(network_data["bus"])

inflation_factors = [0.05, 0.1, 0.2, 0.3, 0.4]


vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_real0.05.mat")
vars_jac = matread("../ptdf_matrices.mat")

pp = vars["coeff_p"]
active_buses = [i for i=1:num_bus if norm(pp[:,i])>1e-5]
num_active = length(active_buses)

pp_jac = vars_jac["Hac_f"][1:38,1:24]
pq_jac = vars_jac["Hac_f"][1:38,25:48]
qp_jac = vars_jac["Hac_f"][39:76,1:24]
qq_jac = vars_jac["Hac_f"][39:76,25:48]

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

    # p = [pp[:,active_buses] qp[:,active_buses]]
    # p_jac = [pp_jac[:,active_buses] qp_jac[:,active_buses]]

    # savefig(heatmap(1:38,2*active_buses,abs(p' - p_jac'),aspect_ratio=1),"plots"string(num_bus)"/heatmaps/heatmap_p_"string(inflation)".pdf")
    # savefig(heatmap(1:38,2*active_buses,,aspect_ratio=1),"heatmap_q_"string(inflation)".pdf")


    # savefig(heatmap(1:num_active,1:38,abs(pp[:,active_buses]-pp_jac[:,active_buses]),aspect_ratio=1, xaxis=(0:4:20)),"plots"string(num_bus)"/heatmaps/heatmap_pp_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:num_active,1:38,abs(pq[:,active_buses]-pq_jac[:,active_buses]),aspect_ratio=1, xaxis=(0:4:20)),"plots"string(num_bus)"/heatmaps/heatmap_pq_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:num_active,1:38,abs(qp[:,active_buses]-qp_jac[:,active_buses]),aspect_ratio=1, xaxis=(0:4:20)),"plots"string(num_bus)"/heatmaps/heatmap_qp_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:num_active,1:38,abs(qq[:,active_buses]-qq_jac[:,active_buses]),aspect_ratio=1, xaxis=(0:4:20)),"plots"string(num_bus)"/heatmaps/heatmap_qq_"string(convert(Int64,inflation*100))".pdf")

    savefig(heatmap(1:38,1:num_active,abs(pp[:,active_buses]-pp_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_pp_"string(convert(Int64,inflation*100))".pdf")
    savefig(heatmap(1:38,1:num_active,abs(pq[:,active_buses]-pq_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_pq_"string(convert(Int64,inflation*100))".pdf")
    savefig(heatmap(1:38,1:num_active,abs(qp[:,active_buses]-qp_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_qp_"string(convert(Int64,inflation*100))".pdf")
    savefig(heatmap(1:38,1:num_active,abs(qq[:,active_buses]-qq_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_qq_"string(convert(Int64,inflation*100))".pdf")
end

# for i = 1:num_branch
#     savefig(plot(inflation_factors,p_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_real_line_"string(i)".pdf")
#     savefig(plot(inflation_factors,q_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_reactive_line_"string(i)".pdf")
# end
