include("./integration_meth.jl")
using Plots

euler = MathConstants.e
pi = MathConstants.pi

bs = [0, 0.5, 0.5, 1]
cs = [1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]
as = [0.0  0.0  0.0   0.0;
      0.5  0.0  0.0   0.0;
      0.0  0.5  0.0   0.0;
      0.0  0.0  1.0   0.0]

# start values
t_start = 0
t_end = 5
h = 0.1

x_start = pi/2
dx_start = p0 = 0
X_start = [x_start, dx_start]

# analytical solution
x(t) = 1/144*t.^4 .+ pi/2

# system of first order differential equations
# https://www.youtube.com/watch?v=_xA-O3TazWQ
f(x) = sqrt(x - x_start)
p(x) = sqrt(4/3 * (x - x_start)^(3/2))

dp(x) = f(x)

P(X, t) = [p(X[1]), f(X[1])]


# evaluate
ts = collect(t_start:h:t_end)
xs_analy = x(ts)
xs_euler = expl_euler_vec_all(t_start, t_end, X_start, P, h)
xs_rk = standard_runge_kutta_vec_all_foo(t_start, t_end, X_start, P, h)



# alternative start value
t_start_alt = 0.01
x_start_alt = x(t_start_alt)
dx_start_alt = p(x_start_alt)
X_start_alt = [x_start_alt, dx_start_alt]

ts_alt = collect(t_start_alt:h:t_end)
xs_euler_alt = expl_euler_vec_all(t_start_alt, t_end, X_start_alt, P, h)
xs_rk_alt = standard_runge_kutta_vec_all(t_start_alt, t_end, X_start_alt, P, h)



#plot
plotly()
pl = plot(ts, xs_analy, title="x(t) = 1/144*t^4 + pi/2, h="*string(h), label="analytical", xlabel="t", ylabel="x(t)", size=(800,600))
plot!(ts, xs_euler[:,1], label="explicit euler")
plot!(ts, xs_rk[:,1], label="runge kutta")
plot!(ts_alt, xs_euler_alt[:,1], label="explicit euler t0="*string(t_start_alt))
plot!(ts_alt, xs_rk_alt[:,1], label="runge kutta t0="*string(t_start_alt))
savefig(pl, "../plot/e6_t4.html")
#display(p)

pyplot()
pl = plot(ts, xs_analy, title="x(t) = 1/144*t^4 + pi/2, h="*string(h), label="analytical", xlabel="t", ylabel="x(t)", size=(800,600))
plot!(ts, xs_euler[:,1], label="explicit euler")
plot!(ts, xs_rk[:,1], label="runge kutta")
plot!(ts_alt, xs_euler_alt[:,1], label="explicit euler t0="*string(t_start_alt))
plot!(ts_alt, xs_rk_alt[:,1], label="runge kutta t0="*string(t_start_alt))
savefig(pl, "../plot/e6_t4.pdf")
