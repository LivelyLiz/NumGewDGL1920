#Task 3
include("./integration_meth.jl")
using LinearAlgebra
using Plots
pyplot()

euler = MathConstants.e

#start values
h = 0.001
x_start = 0
t_start = 0
t_end = 10.0

#functions
f(x, t) = -x + euler^(-t)*cos(t)
df(x, t) = -1.0
x(t) = euler.^(-t).*sin.(t)

function impl_linear_two_step_all(t_start, t_end, x_start, f, df, h, alpha, newton_steps=1)
    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1:2] = standard_runge_kutta_all(t_start, t_start+h, x_start, f, h)

    alphas = [1.0, -alpha-1.0, alpha] #a_2, a_1, a
    betas = [-5.0*alpha-1.0, -8.0*alpha+8.0, 1.0*alpha+5.0].*1.0/12.0 #b_0, b_1, b_2

    for i = 2:(size(ts,1)-1)
        #implicit two step method
        #x_j+1 = -a_1*x_j -a*x_j-1 + h*(b_2*f_j+1 + b_1*f_j + b_0*f_j-1)
        x_approx = x_res[i]
        x_a = alphas[2]*x_res[i] + alphas[3]*x_res[i-1]
        x_b = h*(betas[2]*f(x_res[i], ts[i]) + betas[1]*f(x_res[i-1], ts[i-1]))
        for j = 1:newton_steps
            #use newton's method to approximate 0 = x_j+1 - x_j - h*f(x_j+1, t_j+1)
            #newton: a_n+1 = a_n - f(a_n)/f'(a_n)
            x_approx = x_approx - (alphas[1]*x_approx-h*betas[3]*f(x_approx,ts[i+1])+x_a-x_b)/(1-h*betas[3]*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_approx
    end

    return x_res
end

alphas = [-1.0, -0.999, -0.9, -0.8, -0.3, 0, 0.3, 0.8, 0.9, 0.999, 1.0]

ts = collect(t_start:h:t_end)
x_analy = x(ts)
p = plot(ts, x_analy, label = "Analytical", width=2, linealpha=0.5, size=(1000, 800))

for alpha in alphas
    xs_impl = impl_linear_two_step_all(t_start, t_end, x_start, f, df, h, alpha)
    println("Error for alpha $alpha, h $h, t_end $t_end:\n $(xs_impl[end]-x_analy[end])")
    plot!(p, ts, xs_impl, label = "alpha = $alpha")
end

savefig(p, "../plot/e10_t3.html")
display(p)
