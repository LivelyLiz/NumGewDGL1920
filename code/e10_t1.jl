#Task 1
include("./integration_meth.jl")
using LinearAlgebra

euler = MathConstants.e

#start values
h = 10.0^(-3.0)
X_start = [1; 0; 0]
t_start = 0
t_end = 0.3

X_analy = [9.886739393819 * 10.0^(-1.0);
           3.447715743689 * 10.0^(-5.0);
           1.129158346063 * 10.0^(-2.0)]

#functions
f(X, t) = [-0.04*X[1]+10.0^4.0*X[2]*X[3];
            0.04*X[1]-10.0^4.0*X[2]*X[3]-3*10.0^7.0*X[2]^2.0;
                                         3*10.0^7.0*X[2]^2.0]
df(X, t) = [-0.04   10.0^4.0*X[3]                     10.0^4.0*X[2];
             0.04  -10.0^4.0*X[3]-6*10.0^7.0*X[2]    -10.0^4.0*X[2];
             0.0    6*10.0^7.0*X[2]                   0.0]

newton_steps = [1, 10, 20]

# first calls to time need so compile the functions so this is time warmup
# https://docs.julialang.org/en/latest/manual/performance-tips/#Measure-performance-with-[@time](@ref)-and-pay-attention-to-memory-allocation-1
println("Warmup calls")
@time standard_runge_kutta_vec_all(t_start, t_end, X_start, f, h)
@time impl_euler_vec_all(t_start, t_end, X_start, f, df, h, 1)

# explicit single step
println("Runge Kutta:\nTime:")
xs_ess = @time standard_runge_kutta_vec_all(t_start, t_end, X_start, f, h)
println("Error: \n $(xs_ess[end,:]-X_analy)\n")

# implicit euler with newton
for i in newton_steps
    println("Implicit Euler with $i Newton steps:\nTime:")
    xs_ie = @time impl_euler_vec_all(t_start, t_end, X_start, f, df, h, i)
    println("Error: $(xs_ie[end,:]-X_analy)\n")
end
