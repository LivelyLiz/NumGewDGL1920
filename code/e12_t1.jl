include("integration_meth.jl")

using Plots
gr()

# start values
t_start = 0.0
t_end = 10.0
y_start = [1; 1]
h_start = 0.01

# adaptive rk values
tolerance = 10^(-4)
as = [0 0;
      1 0]
bs = [0.5 0.5;
      1.0 0.0]
cs = [0, 1]

# function
f(y, t) = [y[2] -10*(y[1]^2-1)*y[2]-y[1]]

# result
ys, ts = adaptive_expl_runge_kutta_vec_all(t_start, t_end, y_start, f, h_start, 2, bs, cs, as, tolerance, 1)

solution_plot = plot(ts, [ys[i][1] for i in 1:length(ys)], label="Van-der-Pol", xlabel="t", ylabel="x")
h_plot = plot(ts[2:length(ts)], [ts[i] - ts[i-1] for i in 2:length(ts)], label="h_i", xlabel = "t", ylabel = "h")
p = plot(solution_plot, h_plot)
savefig(p, "../plot/e12_t1.pdf")
