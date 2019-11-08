#Task 2
include("./integration_meth.jl")
using LinearAlgebra

euler = MathConstants.e

#start values
h = 0.1
x_start = [1, 0]
t_start = 0
t_end = 1

#functions
dX(X,t) = [-X[1]+X[2], X[1]-X[2]]
X(t) = [euler^(-2t)/2 + 0.5, -euler^(-2t)/2 + 0.5]


ts = collect(t_start:h:t_end)
real_vals = [X(t) for t in t_start:h:t_end]
expl_euler_vals = expl_euler_vec_all(t_start, t_end, x_start, dX, h)

L = 2.0
error_euklid = 0.0

for i in 1:Int(((t_end-t_start)/h))
    upper = i * h
    lower = (i-1) * h

    c = norm(dX(expl_euler_vals[i,:], i))
    integral = -1/2*euler^((t_end-upper)*L)*h - 1/4 * euler^((t_end-upper)*L) + 1/4 * euler^((t_end-lower)*L)

    global error_euklid += c*integral
end

error_log = 0

for i in 1:Int(((t_end-t_start)/h))
    upper = i * h
    lower = (i-1) * h

    c = norm(dX(expl_euler_vals[i,:], i))
    integral = 1/2 * h^2

    global error_log += c*integral
end

println("Explicit euler values")
display(expl_euler_vals)
print("\nError with euclidean norm   : "*string(error_euklid))
print("\nError with logarithmic norm : "*string(error_log))
