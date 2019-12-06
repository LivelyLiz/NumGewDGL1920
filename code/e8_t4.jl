include("./integration_meth.jl")
using Plots
using LinearAlgebra

euler = MathConstants.e

# start values
t_start = 0
t_end_E = 31557600
t_end_M = 59355072

n = 200000.0

r_E0  = [150.0*10.0^9; 0]
dr_E0 = [0;          29.0 * 10.0^3]

r_M0  = [228.0*10.0^9; 0]
dr_M0 = [0;          24.0 * 10.0^3]

r_S0  = [0; 0]
dr_S0 = [0; 0]

X_start = [ r_E0[1];  #1
            r_E0[2];  #2
           dr_E0[1];  #3
           dr_E0[2];  #4
            r_M0[1];  #5
            r_M0[2];  #6
           dr_M0[1];  #7
           dr_M0[2];  #8
            r_S0[1];  #9
            r_S0[2];  #10
           dr_S0[1];  #11
           dr_S0[2];] #12

#constants
gamma = 6.672 * 10.0^(-11)
m_E = 5.98 * 10.0^24
m_M = 6.42 * 10.0^23
m_S = 1.99 * 10.0^30

print(m_E)

gem = gamma * m_E * m_M
ges = gamma * m_E * m_S
gms = gamma * m_M * m_S

#functions
r_E(X) = [X[1]; X[2]]
r_M(X) = [X[5]; X[6]]
r_S(X) = [X[9]; X[10]]

part(r_a, r_b) = (r_a-r_b)/(norm(r_a-r_b))^3

F(X, t) = [ X[3];
            X[4];
           (-gem * part(r_E(X), r_M(X))[1] + ges * part(r_S(X), r_E(X))[1])/m_E;
           (-gem * part(r_E(X), r_M(X))[2] + ges * part(r_S(X), r_E(X))[2])/m_E;
            X[7];
            X[8];
            (gem * part(r_E(X), r_M(X))[1] + gms * part(r_S(X), r_M(X))[1])/m_M;
            (gem * part(r_E(X), r_M(X))[2] + gms * part(r_S(X), r_M(X))[2])/m_M;
            X[11];
            X[12];
           (-ges * part(r_S(X), r_E(X))[1] - gms * part(r_S(X), r_M(X))[1])/m_S;
           (-ges * part(r_S(X), r_E(X))[2] - gms * part(r_S(X), r_M(X))[2])/m_S
           ]

#evaluation

xs_euler = expl_euler_vec_all(t_start, t_end_M, X_start, F, t_end_M/n)
xs_erk = standard_runge_kutta_vec_all(t_start, t_end_M, X_start, F, t_end_M/n)
xs_ab = adams_bashforth_5_vec_all(t_start, t_end_M, X_start, F, t_end_M/n)

#plot
plotly()
p = plot(xs_euler[:,1], xs_euler[:,2],label="Euler Earth", xlabel="x(t)", ylabel="y(t)")
plot!(p, xs_erk[:,1], xs_erk[:,2],label="ERK Earth")
plot!(p, xs_ab[:,1], xs_ab[:,2],label="Adam Bashforth Earth")
plot!(p, xs_euler[:,5], xs_euler[:,6],label="Euler Mars")
plot!(p, xs_erk[:,5], xs_erk[:,6],label="ERK Mars")
plot!(p, xs_ab[:,5], xs_ab[:,6],label="Adam Bashforth Mars")
plot!(p, xs_euler[:,9], xs_euler[:,10],label="Euler Sun")
plot!(p, xs_erk[:,9], xs_erk[:,10],label="ERK Sun")
plot!(p, xs_ab[:,9], xs_ab[:,10],label="Adam Bashforth Sun")
savefig(p, "../plot/e8_t4.html")

pyplot()
p = plot(xs_euler[:,1], xs_euler[:,2],label="Euler Earth", xlabel="x(t)", ylabel="y(t)")
plot!(p, xs_erk[:,1], xs_erk[:,2],label="ERK Earth")
plot!(p, xs_ab[:,1], xs_ab[:,2],label="Adam Bashforth Earth")
plot!(p, xs_euler[:,5], xs_euler[:,6],label="Euler Mars")
plot!(p, xs_erk[:,5], xs_erk[:,6],label="ERK Mars")
plot!(p, xs_ab[:,5], xs_ab[:,6],label="Adam Bashforth Mars")
plot!(p, xs_euler[:,9], xs_euler[:,10],label="Euler Sun")
plot!(p, xs_erk[:,9], xs_erk[:,10],label="ERK Sun")
plot!(p, xs_ab[:,9], xs_ab[:,10],label="Adam Bashforth Sun")
savefig(p, "../plot/e8_t4.pdf")
