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

inflation_factors = [0.05, 0.1, 0.2, 0.3, 0.4]
algorithms = [1,2]

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



for algo in algorithms
    ctr = 0
    for inflation in inflation_factors
        ctr = ctr+1
        vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_real"string(inflation)"_algorithm_"string(algo)".mat")
        pp = vars["coeff_p"]
        pq = vars["coeff_q"]
        p_err = vars["approx_error"]
        p_error[:,ctr] = p_err
        vars = matread("../results"string(num_bus)"/matrix_forms/linear_approximations_reactive"string(inflation)"_algorithm_"string(algo)".mat")
        qp = vars["coeff_p"]
        qq = vars["coeff_q"]
        q_err = vars["approx_error"]
        q_error[:,ctr] = q_err
    end



    for i = 1:num_branch
        realfig = plot(inflation_factors,p_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",label="Algorithm "string(algo))
        # plot!(inflation_factors,p_error_old[i,:],color=:red,label="gd")
        savefig(realfig,"plots"string(num_bus)"/approximation_error/approx_error_real_line_"string(i)"_algorithm_"string(algo)".pdf")

        reactivefig = plot(inflation_factors,q_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",label="Algorithm "string(algo))
        # plot!(inflation_factors,q_error_old[i,:],color=:red,label="gd")
        savefig(realfig,"plots"string(num_bus)"/approximation_error/approx_error_reactive_line_"string(i)"_algorithm_"string(algo)".pdf")
    end

end
