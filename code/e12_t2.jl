include("integration_meth.jl")

using Plots
pyplot()

# start values
t_start = 0.0
t_end = 17.06522
y_starts = [[0.994; 0; 0; -2.001585106],
            [0.994; 0; 0; -2.01],
            [0.994; 0; 0; -2.02],
            [0.994; 0; 0; -2.03],
            [0.994; 0; 0; -2.031732630]]
h_start = 1

# adaptive rk values
tolerance = 10^(-6)
as = [0 0 0 0 0 0 0;
      0.2 0 0 0 0 0 0;
      3.0/40 9.0/40 0 0 0 0 0;
      44.0/45 -56.0/15 32.0/9 0 0 0 0;
      19372.0/6561 -25360.0/2187 64448.0/6561 -212.0/729 0 0 0;
      9017.0/3168 -355.0/33 46732.0/5247 49.0/176 -5103.0/18656 0 0;
      35.0/384 0 500.0/1113 125.0/192 -2187.0/6784 11.0/84 0]
bs = [35.0/384 0 500.0/1113 125.0/192 -2187.0/6784 11.0/84 0;
      5179.0/57600 0.0 7571.0/16695 393.0/640 -92097.0/339200 187.0/2100 1.0/40]
cs = [0, 0.2, 0.3, 0.8, 8.0/9, 1, 1]

# function
mu = 0.012277471
N1(y) = ((y[1]+mu)^2.0 + y[3]^2.0)^1.5
N2(y) = ((y[1]-(1-mu))^2.0 + y[3]^2.0)^1.5
ddx(a, b, c, d, mu, n1, n2) = a + 2b - (1-mu)*c/n1 - mu*d/n2
f_mu(y, t, mu) = [y[2] ddx(y[1], y[4], y[1]+mu, y[1]-(1-mu), mu, N1(y), N2(y)) y[4] ddx(y[3], -y[2], y[3], y[3], mu, N1(y), N2(y))]
f(y, t) = f_mu(y, t, mu)

# test time
println("Warm-Up")
@time adaptive_expl_runge_kutta_vec_all(t_start, 0.1, y_starts[1], f, h_start, 7, bs, cs, as, tolerance, 4, h_start/1000)
@time standard_runge_kutta_vec_all(t_start, 0.1, y_starts[1], f, 0.001)

colors = ["red", "blue", "green", "magenta", "orange"]

p_s_ad = plot(title="Result adaptive RK")
p_h_ad = plot(title="h Plot")
p_s_rk4 = plot(title="Result RK4")
# result
for (i, y_start) in enumerate(y_starts)
    println("Adaptive RK $(y_start[4])")
    ys, ts = @time adaptive_expl_runge_kutta_vec_all(t_start, t_end, y_start, f, h_start, 7, bs, cs, as, tolerance, 4, false)
    println("Standard RK $(y_start[4])")
    ys_rk4 = @time standard_runge_kutta_vec_all(t_start, t_end, y_start, f, 0.001)
    plot!(p_s_ad, [ys[i][1] for i in 1:length(ys)], [ys[i][3] for i in 1:length(ys)], label="$(y_start[4])", color=colors[i], aspect_ratio=1)
    plot!(p_h_ad, ts[2:length(ts)], [ts[i] - ts[i-1] for i in 2:length(ts)], label="h_i", xlabel = "t", ylabel = "h", color=colors[i])
    plot!(p_s_rk4, ys_rk4[:,1], ys_rk4[:,3], label="$(y_start[4])", color=colors[i], aspect_ratio=1)
end
p = plot(p_s_ad, p_h_ad, p_s_rk4, size=(1024, 1024))
savefig(p, "../plot/e12_t2.pdf")
