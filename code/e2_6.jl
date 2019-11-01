#Task 6
include("./integration_meth.jl")
using Plots

euler = MathConstants.e

#start values
k_max = 12
t_start = 0
t_end = 1/2
x_start = 1

#functions
x(t) = -1.0/(t-1)
f(x, t) = x^2
err_est(h) = 4*(euler^2-1)*h

#init
approx = zeros(k_max)
estErr = zeros(k_max)
numErr = zeros(k_max)

#compute for ks
for k = 1:k_max
    h = 2.0^(-k)
    approx[k] = expl_euler(t_start, t_end, x_start, f, h) #approximated value at x(1/2)
    estErr[k] = err_est(h)                                #estimated error
    numErr[k] = abs(approx[k]-2)                          #numerical error
end

#plot
ks = collect(1:k_max)

plotly()
p = plot(ks, estErr, label="error estimated", xlabel="k", xticks = 0:1:k_max, ylabel="error")
plot!(p, ks, numErr, label="error numerical")
savefig(p, "../plot/e2_t6.html")
#display(p)

gr()
p = plot(ks, estErr, label="error estimated", xlabel="k", xticks = 0:1:k_max, ylabel="error")
plot!(p, ks, numErr, label="error numerical")
savefig(p, "../plot/e2_t6.pdf")
