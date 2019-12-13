#Task 2
include("./integration_meth.jl")
#using Plots
using LinearAlgebra
using Plots
plotly()

euler = MathConstants.e

#start values
hs = [0.25, 0.2, 0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.0025, 0.002, 0.001, 0.0005, 0.00025, 0.0002, 0.0001]
x_start = euler
t_start = 0
t_end = 1

#functions
f(x, t) = -t^2.0 * x
df(x, t) = -t^2.0
x(t) = euler.^(1 .- (t.^3)./3)

# inital value methods
init_val_meths = [expl_euler_all,
                   second_order_runge_kutta_all,
                   heun_runge_kutta_all,
                   standard_runge_kutta_all]

init_val_meth_names = ["Explicit Euler",
                       "2nd order RK",
                       "Heun",
                       "Stand. RK"]

# integration methods
int_meth_names = ["Adams-Bashf.",
                  "Adams-Moulton",
                  "Nyström",
                  "Milne-Simpson"]

function getVals_integration_method(method, init_val_meth, t_start, t_end, x_start, f, df, h)

    if method == 1
        return adams_bashforth_all(t_start, t_end, x_start, f, h, 3, init_val_meth)
    elseif method == 2
        return adams_moulton_all(t_start, t_end, x_start, f, df, h, 3, init_val_meth)
    elseif method == 3
        return nystroem_all(t_start, t_end, x_start, f, h, 3, init_val_meth)
    elseif method == 4
        return milne_simpson_all(t_start, t_end, x_start, f, df, h, 2, init_val_meth)
    else
        return
    end

end

# wanted order
orders = [3, 4, 3, 4]

ts = collect(t_start:hs[end]:t_end)
xs_analy = x(ts)

subplots = []

for int_meth in 1:4
    subplot = plot(title=int_meth_names[int_meth], xlabel = "index of h", ylabel="ooc", xticks=1:1:length(hs)-1)
    plot!(subplot, [1, length(hs)-1], [orders[int_meth], orders[int_meth]], label="expected order", width=2, linealpha=0.5, linestyle=:dot)

    for init_val_meth in 1:4
        prev_sol = []
        sol = []
        or_o_con = []
        for (i,h) in enumerate(hs)
            if i == 1
                sol = getVals_integration_method(int_meth, init_val_meths[init_val_meth], t_start, t_end, x_start, f, df, h)
                continue
            end

            prev_sol = sol
            sol = getVals_integration_method(int_meth, init_val_meths[init_val_meth], t_start, t_end, x_start, f, df, h)

            err_1 = abs(prev_sol[end] - xs_analy[end])
            err_2 = abs(sol[end] - xs_analy[end])

            push!(or_o_con, ooc(hs[i-1], h, err_1, err_2))
        end
        plot!(subplot, collect(1:length(hs)-1), or_o_con, label=init_val_meth_names[init_val_meth])

        println(int_meth_names[int_meth]*" "*init_val_meth_names[init_val_meth])
        println(or_o_con)
    end
    push!(subplots, subplot)
end

main_plot = plot(subplots..., layout=(2,2), size=(1024, 800))
savefig("../plot/e9_t1.html")
display(main_plot)
