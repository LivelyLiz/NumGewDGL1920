#Task 2
include("./integration_meth.jl")
include("./e5_t4.jl")
using Plots

euler = MathConstants.e

#constants
sigma = 0.07274
g = 9.81
roh = 998.2
z0 = -1 * 10^(-3)

#start values
s_start = 0
ss_end = collect(0.001:0.001:0.005)
r_start = 0
z_start = z0
phi_start = 0
x_start = [r_start, z_start, phi_start]

#functions
dr(x, s) = cos(x[3])
dz(x, s) = sin(x[3])
function dphi_ks(x, s, ks)
    if (x[1] > 0)
        return 2*ks - (roh*g)/sigma*x[2] - sin(x[3])/x[1]
    else (x[1] <= 0)
        return ks - (roh*g)/(2*sigma)*x[2]
    end
end

F_ks(x, s, ks) = [dr(x,s), dz(x,s), dphi_ks(x,s,ks)]

#Task 2
#init
kss = collect(200:100:500)
n = 100

plotly()
send_plots = []

for (j, s_end) in enumerate(ss_end)
    h = (s_end-s_start)/n
    pl = plot(title="s_end = "*string(s_end))
    for (i, ks) in enumerate(kss)
        F(x, s) = F_ks(x, s, ks)
        vals = expl_euler_vec_all(s_start, s_end, x_start, F, h)
        plot!(pl, vals[:,1], vals[:,2], label="ks = "*string(ks))
    end
    push!(send_plots, pl)
end

plot_euler = plot(send_plots..., size=(1024,1024), xlabel="r(s) in m", ylabel="z(s) in m")
savefig(plot_euler, "../plot/e5_t2_expl_euler.html")

#Task 4
#init
h = 0.1 * 10^(-3) #because we use meters instead of mm

send_plots = []

for (j, s_end) in enumerate(ss_end)
    h = (s_end-s_start)/n
    pl = plot(title="s_end = "*string(s_end))
    for (i, ks) in enumerate(kss)
        F(x, s) = F_ks(x, s, ks)
        vals = runge_kutta_vec_all(s_start, s_end, x_start, F, h)
        plot!(pl, vals[:,1], vals[:,2], label="ks = "*string(ks))
    end
    push!(send_plots, pl)
end

plot_euler = plot(send_plots..., size=(1024,1024), xlabel="r(s) in m", ylabel="z(s) in m")
savefig(plot_euler, "../plot/e5_t4_runge_kutta.html")
