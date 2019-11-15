#Task 1
include("./integration_meth.jl")
using Plots

euler = MathConstants.e

#start values
x_start = 1.0
x_end = 2.0
y_start = 1.0/2.0
h = 0.1

#functions
y(x) = -x./2 .- 1/4 .+ 5/4*euler.^(2*x .- 2)
f(y, x) = x + 2.0*y
f2(dy, x) = 1.0 + 2.0*dy

#order of convergence
function ooc(h1, h2, err1, err2)
    return (log(err1)-log(err2))/(log(h1)-log(h2))
end

#init
xs_approx = collect(x_start:h:x_end)
xs_analy = collect(x_start:0.00001:x_end)

ys_analy = y(xs_analy)
ys_approx = expl_euler_all(x_start, x_end, y_start, f, h)
ys_approx_mod = expl_euler_mod_all(x_start, x_end, y_start, f, f2, h)

#numerical order of convergence
ks = collect(1:5)
errs = zeros(size(ks,1), 2)

for i in 1:size(ks, 1)
    h_i = 10.0^(-ks[i])
    val_expl = expl_euler(x_start, x_end, y_start, f, h_i)
    val_mod = expl_euler_mod(x_start, x_end, y_start, f, f2, h_i)

    err_expl = abs(ys_analy[end] - val_expl)
    err_mod = abs(ys_analy[end] - val_mod)

    errs[i, 1] = err_expl
    errs[i, 2] = err_mod

    if(i > 1)
        p_expl = ooc(10.0^(-ks[i-1]), h_i, errs[i-1,1], errs[i, 1])
        p_mod = ooc(10.0^(-ks[i-1]), h_i, errs[i-1,2], errs[i, 2])

        println("h_i: "*string(h_i))
        println("conv expl: "*string(p_expl))
        println("conv mod: "*string(p_mod))
    end
end

#plot
plotly()
p = plot(xs_approx, ys_approx, label="expl. euler", xlabel="x", xticks = x_start:h:x_end, ylabel="y(x)", yticks=0:0.5:8, size=(800,600))
plot!(xs_approx, ys_approx_mod, label="expl. euler mod")
plot!(xs_analy, ys_analy, label="analytical solution", width=2)
savefig(p, "../plot/e5_t1.html")
#display(p)

gr()
p = plot(xs_approx, ys_approx, label="expl. euler", xlabel="x", xticks = x_start:h:x_end, ylabel="y(x)", yticks=0:0.5:8, size=(800,600))
plot!(xs_approx, ys_approx_mod, label="expl. euler mod")
plot!(xs_analy, ys_analy, label="analytical solution", width=2)
savefig(p, "../plot/e5_t1.pdf")
