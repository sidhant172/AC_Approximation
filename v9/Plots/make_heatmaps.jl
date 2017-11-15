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

inflation_factors = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]


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

pp_dict = Dict{String,Any}()
pq_dict = Dict{String,Any}()
qp_dict = Dict{String,Any}()
qq_dict = Dict{String,Any}()
const_dict = Dict{String,Any}()


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

    pp_dict[string(convert(Int64,inflation*100))] = pp
    pq_dict[string(convert(Int64,inflation*100))] = pq
    qp_dict[string(convert(Int64,inflation*100))] = qp
    qq_dict[string(convert(Int64,inflation*100))] = qq

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

    # savefig(heatmap(1:38,1:num_active,abs(pp[:,active_buses]-pp_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_pp_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:38,1:num_active,abs(pq[:,active_buses]-pq_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_pq_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:38,1:num_active,abs(qp[:,active_buses]-qp_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_qp_"string(convert(Int64,inflation*100))".pdf")
    # savefig(heatmap(1:38,1:num_active,abs(qq[:,active_buses]-qq_jac[:,active_buses])',aspect_ratio=1, yaxis=(0:2:20)),"plots"string(num_bus)"/heatmaps/heatmap_qq_"string(convert(Int64,inflation*100))".pdf")

    ###### making rms maps instead



end

# combined_rms_plot = plot(size = (1400,800), layout=4, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, linewidth=2,
#     lab = "Active power lp", title_position = :bottom, xtickfont = font(8, "Courier"), xlabel = "Radius", ylabel = "RMS(Optimal-Taylor)", inflation_factors, (1/sqrt(38*length(active_buses)))*[norm(pp_dict[string(convert(Int64,100*inflation))][:,active_buses]-pp_jac[:,active_buses]) for inflation in inflation_factors]  )
# plot!(right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, linewidth=2, lab = "Active power lq", xtickfont = font(8, "Courier"), xlabel = "Radius", ylabel = "RMS(Optimal-Taylor)", inflation_factors, (1/sqrt(38*length(active_buses)))*[norm(pq_dict[string(convert(Int64,100*inflation))][:,active_buses]-pq_jac[:,active_buses]) for inflation in inflation_factors], subplot = 2 )
# plot!(right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, linewidth=2, color = :red, lab = "Reactive power lp", xtickfont = font(8, "Courier"), xlabel = "Radius", ylabel = "RMS(Optimal-Taylor)", inflation_factors, (1/sqrt(38*length(active_buses)))*[norm(qp_dict[string(convert(Int64,100*inflation))][:,active_buses]-qp_jac[:,active_buses]) for inflation in inflation_factors], subplot = 3 )
# plot!(right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, linewidth=2, color = :red, lab = "Reactive power lq", xtickfont = font(8, "Courier"), xlabel = "Radius", ylabel = "RMS(Optimal-Taylor)", inflation_factors, (1/sqrt(38*length(active_buses)))*[norm(qq_dict[string(convert(Int64,100*inflation))][:,active_buses]-qq_jac[:,active_buses]) for inflation in inflation_factors], subplot = 4  )
# savefig(combined_rms_plot,"combined_rms_plot.pdf")


combined_rms_plot = plot(size = (1400,350), layout=2, right_margin=40px, left_margin=20px, top_margin= 10px, bottom_margin=10px, linewidth=2,
    lab = "Real power lp", title_position = :bottom, xtickfont = font(8, "Courier"), xlabel = "Radius (percentage)", ylabel = "RMS(Optimal-Taylor)", convert(Array{Int64,1},100*inflation_factors), (1/sqrt(38*length(active_buses)))*[norm(pp_dict[string(convert(Int64,100*inflation))][:,active_buses]-pp_jac[:,active_buses]) for inflation in inflation_factors]  )
plot!(right_margin=40px, left_margin=20px, top_margin= 10px, bottom_margin=10px, linewidth=2, linestyle=:dash, color = :red, lab = "Real power lq", xtickfont = font(8, "Courier"), xlabel = "Radius (percentage)", ylabel = "RMS(Optimal-Taylor)", convert(Array{Int64,1},100*inflation_factors), (1/sqrt(38*length(active_buses)))*[norm(pq_dict[string(convert(Int64,100*inflation))][:,active_buses]-pq_jac[:,active_buses]) for inflation in inflation_factors], subplot = 1 )
plot!(right_margin=40px, left_margin=20px, top_margin= 10px, bottom_margin=10px, linewidth=2, lab = "Reactive power lp", xtickfont = font(8, "Courier"), xlabel = "Radius (percentage)", ylabel = "RMS(Optimal-Taylor)", convert(Array{Int64,1},100*inflation_factors), (1/sqrt(38*length(active_buses)))*[norm(qp_dict[string(convert(Int64,100*inflation))][:,active_buses]-qp_jac[:,active_buses]) for inflation in inflation_factors], subplot = 2 )
plot!(right_margin=40px, left_margin=20px, top_margin= 10px, bottom_margin=10px, linewidth=2, linestyle=:dash, color = :red, lab = "Reactive power lq", xtickfont = font(8, "Courier"), xlabel = "Radius (percentage)", ylabel = "RMS(Optimal-Taylor)", convert(Array{Int64,1},100*inflation_factors), (1/sqrt(38*length(active_buses)))*[norm(qq_dict[string(convert(Int64,100*inflation))][:,active_buses]-qq_jac[:,active_buses]) for inflation in inflation_factors], subplot = 2  )

savefig(combined_rms_plot,"combined_rms_plot.pdf")

# combined_heatmap = heatmap(size = (1400,800), layout=4, right_margin=10px, left_margin=10px, top_margin= 0px, bottom_margin=0px,
#     1:38,1:num_active, abs(pp_dict["5"][:,active_buses]-pp_jac[:,active_buses])'
#     ,aspect_ratio=1, yaxis=(0:2:20))
# heatmap!(1:38,1:num_active, abs(pq_dict["5"][:,active_buses]-pq_jac[:,active_buses])', aspect_ratio=1, yaxis=(0:2:20), subplot=2)
# heatmap!(1:38,1:num_active, abs(pp_dict["40"][:,active_buses]-pp_jac[:,active_buses])', aspect_ratio=1, yaxis=(0:2:20), subplot=3)
# heatmap!(1:38,1:num_active, abs(pq_dict["40"][:,active_buses]-pq_jac[:,active_buses])', aspect_ratio=1, yaxis=(0:2:20), subplot=4)
#
# annotate!(2,18,text("Lp R=0.05",:white,:left), subplot=1)
# annotate!(2,18,text("Lq R=0.05",:white,:left), subplot=2)
# annotate!(2,18,text("Lp R=0.4",:white,:left), subplot=3)
# annotate!(2,18,text("Lq R=0.4",:white,:left), subplot=4)
#
# savefig(combined_heatmap,"combined_heatmap.pdf")






# combined_heatmap = heatmap(size = (1400,900), layout=4, right_margin=20px, left_margin=20px, top_margin= 20px, bottom_margin=10px,
# 1:num_active,1:38,abs(pp_dict["5"][:,active_buses]-pp_jac[:,active_buses]),aspect_ratio=1, "combined_heatmap.pdf" )

# for i = 1:num_branch
#     savefig(plot(inflation_factors,p_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_real_line_"string(i)".pdf")
#     savefig(plot(inflation_factors,q_error[i,:],xlabel="Radius",ylabel="Maximum approximation error",leg=false),"plots"string(num_bus)"/approximation_error/approx_error_reactive_line_"string(i)".pdf")
# end
