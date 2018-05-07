using MAT
using Plots
using Plots.PlotMeasures
gr()


dirname = "../results24"
plot_dirname = "plots24"


vars = matread(string(dirname)"/error_tracking_real_line1.mat")
# vars = matread("../results14/error_tracking_real_line1.mat")
# vars = matread("../results24/error_tracking_real_line1.mat")
# vars = matread("../results57/error_tracking_real_line1.mat")
nlp_err = max(vars["nlp_err_pos"],vars["nlp_err_neg"])
lp_err = vars["lp_err"]

@show ind_phase1 = minimum(find(lp_err.>0))



err_plot = plot(size = (1400,350), layout=2, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px,
    nlp_err[2:end], xlabel = "Iteration number", ylabel="Maximum approximation error", lab="Line 1 Real",  xtickfont = font(10, "Courier"), linewidth=2)

plot!(lp_err, xlabel = "Iteration number", ylabel="LP Objective", xtickfont = font(10, "Courier"),
     linewidth=2, lab="Line 1 Real", subplot=2, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)
# # savefig(err_real,"plots24/error_tracking_real.pdf")


# err_plot = plot(size = (700,350), right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, lab = "Real",  ylims = [0,0.25],
#     ind_phase1+1:length(nlp_err),nlp_err[ind_phase1+1:end], xlabel = "Iteration number", ylabel="Maximum approximation error", xtickfont = font(10, "Courier"), linewidth=2)
#
#
# plot!(lp_err[ind_phase1:end], xlabel = "Iteration number", ylabel="Objective of LP relaxation", color=:red, xtickfont = font(10, "Courier"),
#      linewidth=2, leg=false, subplot=3, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)
# savefig(err_real,"plots24/error_tracking_real.pdf")

lp_deltas_real = vars["lp_deltas"]

# @show string(dirname)"/error_tracking_reactive_line1.mat"
vars = matread(string(dirname)"/error_tracking_reactive_line1.mat")
# vars = matread("../results14/error_tracking_reactive_line1.mat")
# vars = matread("../results24/error_tracking_reactive_line1.mat")
# vars = matread("../results57/error_tracking_reactive_line1.mat")

nlp_err = max(vars["nlp_err_pos"],vars["nlp_err_neg"])
lp_err = vars["lp_err"]

@show ind_phase1 = minimum(find(lp_err.>0))
#
plot!(nlp_err[2:end], xlabel = "Iteration number", linestyle =:dash, ylabel="Maximum approximation error \n (pu)",xtickfont = font(10, "Courier"), linewidth=2, lab="Line 1 Reactive", subplot=1, right_margin=40px, left_margin=20px, top_margin=30px, bottom_margin=10px)
# # plot!(nlp_err)
plot!(lp_err, xlabel = "Iteration number", linestyle = :dash, ylabel="LP Objective (pu)", xtickfont = font(10, "Courier"), linewidth=2, lab="Line 1 Reactive", subplot=2, right_margin=40px, left_margin=20px, top_margin=30px, bottom_margin=10px)
# #
# lp_deltas_reactive = vars["lp_deltas"]

# annotate!(65,1,text("Line #1 Active power", font(9), :left), subplot=1)
# annotate!(95,0.8,text("Line #1 Reactive power", font(9), :left), subplot=2)
# annotate!(65,0.0004,text("Line #1 Active power", font(9), :left), subplot=3)
# annotate!(95,0.0016,text("Line #1 Reactive power", font(9), :left), subplot=4)

# savefig(plot(vars["lp_deltas"]), "plots24/lp_deltas_reactive.pdf")


# plot!(lp_deltas_real, xlabel = "Iteration number", ylabel="Change in lp",xtickfont = font(10, "Courier"), linewidth=2, color=:green, leg=false, subplot=3, right_margin=40, left_margin=20, top_margin= 30, bottom_margin=10)
# plot!(nlp_err)
# plot!(lp_deltas_reactive, xlabel = "Iteration number", ylabel="Change in lq", xtickfont = font(10, "Courier"),color=:green, linewidth=2, leg=false, subplot=4, right_margin=40, left_margin=20, top_margin= 30, bottom_margin=10)


savefig(err_plot, string(plot_dirname)"/error_tracking.pdf")

# savefig(err_plot,"plots14/error_tracking.pdf")
# savefig(err_plot,"plots24/error_tracking.pdf")
# savefig(err_plot,"plots57/error_tracking.pdf")

# plot!(ind_phase1+1:length(nlp_err),nlp_err[ind_phase1+1:end], linestyle = :dash, lab = "Reactive", xlabel = "Iteration number", ylabel="Maximum approximation error",xtickfont = font(10, "Courier") ,linewidth=2, right_margin=40px, left_margin=20px, top_margin=30px, bottom_margin=10px)

# plot!(lp_err[ind_phase1:end], xlabel = "Iteration number", ylabel="Objective of LP relaxation", xtickfont = font(10, "Courier"),color=:red, linewidth=2, leg=false, subplot=4, right_margin=40, left_margin=20, top_margin= 30, bottom_margin=10)


# savefig(err_plot, string(plot_dirname)"/error_tracking_phase1.pdf")



# savefig(plot(lp_deltas_real),"plots24/lp_deltas_real.pdf")
# savefig(plot(lp_deltas_reactive),"plots24/lp_deltas_reactive.pdf")
