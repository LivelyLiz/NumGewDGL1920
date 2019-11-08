#Task 3
include("./integration_meth.jl")
using Plots

euler = MathConstants.e

#start values
hs = [1, 0.1, 0.01, 0.001, 0.00199, 0.002, 0.0021]
y_start = 1
x_start = 0
x_end = 10

#functions
dy(y,x) = -1000*y + 1000*sin(x) + cos(x)
ddy(y,x) = -1000
y(x) = euler^(-1000*x)+sin(x)


#initialization
real_vals = [[y(x) for x in x_start:hs[j]:x_end] for j in 1:length(hs)]
expl_euler_vals = []
impl_euler_vals = []

xs = [collect(x_start:hs[j]:x_end) for j in 1:length(hs)]
plots = []

for (i, h) in enumerate(hs)
    push!(expl_euler_vals, expl_euler_all(x_start, x_end, y_start, dy, h))
    push!(impl_euler_vals, impl_euler_all(x_start, x_end, y_start, dy, ddy, h))
end

#plot
plotly()

for (i, h) in enumerate(hs)
    pl = plot(xs[i], real_vals[i], label="reference", xlabel="x", ylabel="y", title="h="*string(h), ylims=(-10,10))
    plot!(pl, xs[i], expl_euler_vals[i], label="explicit")
    plot!(pl, xs[i], impl_euler_vals[i], label="implicit")
    push!(plots, pl)
end

p = plot(plots..., layout=length(hs), size=(2000,1500))
savefig(p, "../plot/e3_t3.html")
display(p)

pyplot()

plots = []
for (i, h) in enumerate(hs)
    pl = plot(xs[i], real_vals[i], label="reference", title="h="*string(h), ylims=(-10,10))
    plot!(pl, xs[i], expl_euler_vals[i], label="explicit")
    plot!(pl, xs[i], impl_euler_vals[i], label="implicit")
    push!(plots, pl)
end

p = plot(plots..., size=(2000,1500))
savefig(p, "../plot/e3_t3.pdf")
