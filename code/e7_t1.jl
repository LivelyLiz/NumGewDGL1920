include("./integration_meth.jl")
using Plots

euler = MathConstants.e

# start values
t_start_1 = 0
t_end_1 = 4

t_start_2 = 0
t_end_2 = 15

h = 0.01

X_start = [1; 2]

#functions
dX_1(X, t) = [-1   5;
              -5  -1] * X
dX_2(X, t) = [ 0   1;
              -4   0] * X
c1_1 = (1/2-1im)
c2_1 = (1/2+1im)
X_1(t) = [c1_1*euler.^((-1 + 5im).*t)         + c2_1*euler.^((-1 - 5im).*t);
          c1_1*euler.^((-1 + 5im).*t) * 1im   + c2_1*euler.^((-1 - 5im).*t) * -1im]

c1_2 = (1/2-1im/2)
c2_2 = (1/2+1im/2)
X_2(t) = [c1_2*euler.^(2im.*t)                + c2_2*euler.^(-2im.*t);
          c1_2*euler.^(2im.*t) * 2im          + c2_2*euler.^(-2im.*t) * -2im]

#evaluation
ts_1 = collect(t_start_1:h:t_end_1)
ts_2 = collect(t_start_2:h:t_end_2)

Xs_1_analy = X_1(ts_1')
Xs_1_euler = expl_euler_vec_all(t_start_1, t_end_1, X_start, dX_1, h)

Xs_2_analy = X_2(ts_2')
Xs_2_euler = expl_euler_vec_all(t_start_2, t_end_2, X_start, dX_2, h)


#plot
plotly()
p_1 = plot(real.(Xs_1_analy'[:,1]), real.(Xs_1_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
plot!(p_1, Xs_1_euler[:,1], Xs_1_euler[:,2], label="expl. euler")

p_2 = plot(real.(Xs_2_analy'[:,1]), real.(Xs_2_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
plot!(p_2, Xs_2_euler[:,1], Xs_2_euler[:,2], label="expl. euler")
savefig(plot(p_1, p_2), "../plot/e7_t1.html")

pyplot()
p_1 = plot(real.(Xs_1_analy'[:,1]), real.(Xs_1_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
plot!(p_1, Xs_1_euler[:,1], Xs_1_euler[:,2], label="expl. euler")

p_2 = plot(real.(Xs_2_analy'[:,1]), real.(Xs_2_analy'[:,2]), label="analytical (real part)", xlabel="x(t)", ylabel="y(t)", size=(1000,600))
plot!(p_2, Xs_2_euler[:,1], Xs_2_euler[:,2], label="expl. euler")
savefig(plot(p_1, p_2), "../plot/e7_t1.pdf")
