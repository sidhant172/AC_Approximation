using MAT
using Plots
gr()

vars = matread("../results24/error_tracking_real_line1.mat")
# vars = matread("../results57/error_tracking_real_line1.mat")
nlp_err = max(vars["nlp_err_pos"],vars["nlp_err_neg"])
lp_err = vars["lp_err"]



err_plot = plot(size = (1400,800), layout=4, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px, 
    nlp_err[2:end], xlabel = "Iteration number", ylabel="Maximum approximation error", leg=false,  xtickfont = font(10, "Courier"), linewidth=2)

plot!(lp_err, xlabel = "Iteration number", ylabel="Objective of LP relaxation", color=:red, xtickfont = font(10, "Courier"),
     linewidth=2, leg=false, subplot=3, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px,)
# savefig(err_real,"plots24/error_tracking_real.pdf")

lp_deltas_real = vars["lp_deltas"]


vars = matread("../results24/error_tracking_reactive_line1.mat")
# vars = matread("../results57/error_tracking_reactive_line1.mat")

nlp_err = max(vars["nlp_err_pos"],vars["nlp_err_neg"])
lp_err = vars["lp_err"]
plot!(nlp_err[2:end], xlabel = "Iteration number", ylabel="Maximum approximation error",xtickfont = font(10, "Courier"), linewidth=2, leg=false, subplot=2, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)
# plot!(nlp_err)
plot!(lp_err, xlabel = "Iteration number", ylabel="Objective of LP relaxation", xtickfont = font(10, "Courier"),color=:red, linewidth=2, leg=false, subplot=4, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)

lp_deltas_reactive = vars["lp_deltas"]

annotate!(65,1,text("Line #1 Active power", font(9), :left), subplot=1)
annotate!(95,0.8,text("Line #1 Reactive power", font(9), :left), subplot=2)
annotate!(65,0.0004,text("Line #1 Active power", font(9), :left), subplot=3)
annotate!(95,0.0016,text("Line #1 Reactive power", font(9), :left), subplot=4)

# savefig(plot(vars["lp_deltas"]), "plots24/lp_deltas_reactive.pdf")


# plot!(lp_deltas_real, xlabel = "Iteration number", ylabel="Change in lp",xtickfont = font(10, "Courier"), linewidth=2, color=:green, leg=false, subplot=3, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)
# plot!(nlp_err)
# plot!(lp_deltas_reactive, xlabel = "Iteration number", ylabel="Change in lq", xtickfont = font(10, "Courier"),color=:green, linewidth=2, leg=false, subplot=4, right_margin=40px, left_margin=20px, top_margin= 30px, bottom_margin=10px)



savefig(err_plot,"plots24/error_tracking.pdf")
# savefig(err_plot,"plots57/error_tracking.pdf")





# savefig(plot(lp_deltas_real),"plots24/lp_deltas_real.pdf")
# savefig(plot(lp_deltas_reactive),"plots24/lp_deltas_reactive.pdf")
