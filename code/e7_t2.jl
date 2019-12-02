include("./integration_meth.jl")
using Plots

euler = MathConstants.e

# start values
t_start_1 = 0
t_end_1 = 4

t_start_2 = 0
t_end_2 = 15

h_1 = 0.01
h_2 = 0.01

X_start = [1; 2]

#Butcher-Tableaus
names = ["Standard RK", "Methode von England", "3/8 Methode"]
#1 Klassisches RK
#2 Methode von England
as_eng = [0    0    0    0;
          0.5  0    0    0;
          0.25 0.25 0    0;
          0    -1   2    0]
bs_eng = [1/6, 0,   2/3, 1/6]
cs_eng = [0,   1/2, 1/2, 1]

#3 3/8-Regel
as_38 = [0    0    0    0;
          0.5  0    0    0;
          0    1    0    0;
          0    0    1    0]
bs_38 = [1/6, 2/3, 0,   1/6]
cs_38 = [0,   1/2, 1,   1]


#functions
dX_1(X, t) = [-1  5;
              -5  -1] * X
dX_2(X, t) = [ 0  1;
              -4  0] * X

r = 1
c1_1 = (1/2-1im)
c2_1 = (1/2+1im)
X_1(t) = [c1_1*euler.^((-1 + 5im).*t)         + c2_1*euler.^((-1 - 5im).*t);
          c1_1*euler.^((-1 + 5im).*t) * 1im   + c2_1*euler.^((-1 - 5im).*t) * -1im]

c1_2 = (1/2-1im/2)
c2_2 = (1/2+1im/2)
X_2(t) = [c1_2*euler.^(2im.*t)                + c2_2*euler.^(-2im.*t);
          c1_2*euler.^(2im.*t) * 2im          + c2_2*euler.^(-2im.*t) * -2im]

ts_1 = collect(t_start_1:h_1:t_end_1)
ts_2 = collect(t_start_2:h_2:t_end_2)

Xs_1_analy = X_1(ts_1')
Xs_1_rk_stand = standard_runge_kutta_vec_all(t_start_1, t_end_1, X_start, dX_1, h_1)
Xs_1_rk_eng = expl_runge_kutta_vec_all(t_start_1, t_end_1, X_start, dX_1, h_1, 4, bs_eng, cs_eng, as_eng)
Xs_1_rk_38 = expl_runge_kutta_vec_all(t_start_1, t_end_1, X_start, dX_1, h_1, 4, bs_38, cs_38, as_38)
Xs_1 = [Xs_1_rk_stand, Xs_1_rk_eng, Xs_1_rk_38]

Xs_2_analy = X_2(ts_2')
Xs_2_rk_stand = standard_runge_kutta_vec_all(t_start_2, t_end_2, X_start, dX_2, h_2)
Xs_2_rk_eng = expl_runge_kutta_vec_all(t_start_2, t_end_2, X_start, dX_2, h_2, 4, bs_eng, cs_eng, as_eng)
Xs_2_rk_38 = expl_runge_kutta_vec_all(t_start_2, t_end_2, X_start, dX_2, h_2, 4, bs_38, cs_38, as_38)
Xs_2 = [Xs_2_rk_stand, Xs_2_rk_eng, Xs_2_rk_38]

println("Sum Squared Error")
println("1")
sum_sqrd_err_1 = []
for i in 1:size(Xs_1,1)
    push!(sum_sqrd_err_1, sum((Xs_1_analy - Xs_1[i][:,:]').^2, dims=2))
    println(names[i]*" "*string(sum_sqrd_err_1[i]))
end

println("2")
sum_sqrd_err_2 = []
for i in 1:size(Xs_2,1)
    push!(sum_sqrd_err_2, sum((Xs_2_analy - Xs_2[i][:, :]').^2, dims=2))
    println(names[i]*" "*string(sum_sqrd_err_2[i]))
end


#plot
plotly()
p_1 = plot(real.(Xs_1_analy'[:,1]), real.(Xs_1_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
for i in 1:size(names,1)
    plot!(p_1, Xs_1[i][:, 1], Xs_1[i][:, 2], label=names[i])
end

p_2 = plot(real.(Xs_2_analy'[:,1]), real.(Xs_2_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
for i in 1:size(names,1)
    plot!(p_2, Xs_2[i][:, 1], Xs_2[i][:, 2], label=names[i])
end
savefig(plot(p_1, p_2), "../plot/e7_t2.html")

pyplot()
p_1 = plot(real.(Xs_1_analy'[:,1]), real.(Xs_1_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
for i in 1:size(names,1)
    plot!(p_1, Xs_1[i][:, 1], Xs_1[i][:, 2], label=names[i])
end

p_2 = plot(real.(Xs_2_analy'[:,1]), real.(Xs_2_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
for i in 1:size(names,1)
    plot!(p_2, Xs_2[i][:, 1], Xs_2[i][:, 2], label=names[i])
end
savefig(plot(p_1, p_2), "../plot/e7_t2.pdf")
