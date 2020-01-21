"""
X_0   ... start value for newton method
Df    ... Jacobian matrix of multivariate function (Df(X)) taking X as argument
neg_f ... negative of function f taking X as argument
num_iter ... number of iterations for newton method
(t is constant during Newton method)

return the calculated value
"""
function newton_method_vec(X_0, Df, neg_f, num_iter)
    """
    multivariate Newton method
    X_k+1 = X_k - inverse(Df(X_k))*f(X_k)

    to avoid inverse
    solve Df(X_k)*s_k = -f(X_k)
    and X_k+1 = X_k + s_k
    """

    X = X_0
    for i in 1:num_iter
        #backslash operator solving system of equations if possible
        s = Df(X)\neg_f(X)
        X = X + s
    end

    return X
end

"""
Order of convergence
h1...   stepwidth used for error 1
h2...   stepwidth used for error 2
err1... error 1
err2... error 2

Error values are the error between the numerical
and the analytical solution at a certain time t
"""
function ooc(h1, h2, err1, err2)
    return (log(err1)-log(err2))/(log(h1)-log(h2))
end

"""
Estimates the error for adaptive runge kutta method

x_1    ... value (vector) obtained with method of higher Order
xhat_1 ... value obtained with lower order

return highest error over all dimensions
"""
function error_est(x_0, x_1, xhat_1)
    err = 0
    for i in 1:length(x_1)
        val = abs(x_1[i] - xhat_1[i])/abs(1 + x_0[i])
        err = max(err, val)
    end
    return err
end

"""
Calculates the h used for adaptive runge kutta method

tol        ... tolerance
err        ... error (estimated)
min_order  ... minimum order of the methods used
h          ... old value for h

return the new h to use for the ODE solver
"""
function adaptive_h(tol, err, min_order, h)
    return min(2, max(0.5, 0.9*(tol/err)^(1/(min_order+1)))) * h
end

"""
Explicit Euler
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return estimated value of x at point t_end
"""
function expl_euler(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = x_start

    for i = 1:size(ts,1)
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res = x_res + f(x_res, ts[i])*h
    end

    return x_res
end

"""
Explicit Euler with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
f2     ... f2(f(x,t),t) second order derivative of x
h      ... sample point distance

return estimated value of x at point t_end
"""
function expl_euler_mod(t_start, t_end, x_start, f, f2, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = x_start

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j) + f(t_j, f(t_j, x_j))*(t_j+1 - t_j)^2*1/2
        x_res = x_res +
                f(x_res, ts[i])*h +
                f2(f(x_res, ts[i]), ts[i]) * 1/2 * h^2
    end

    return x_res
end

"""
Explicit Euler with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return a vector with inbetween steps
"""
function expl_euler_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros(size(ts,1)+1)
    x_res[1] = x_start

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res[i+1] = x_res[i] + f(x_res[i], ts[i])*h
    end

    return x_res
end

"""
Explicit Euler with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
f2     ... f2(f(x,t),t) second order derivative of x
h      ... sample point distance

return a vector with inbetween steps
"""
function expl_euler_mod_all(t_start, t_end, x_start, f, f2, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros(size(ts,1)+1)
    x_res[1] = x_start

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j) + f(t_j, f(t_j, x_j))*(t_j+1 - t_j)^2*1/2
        #x_res[i+1] = x_res[i] +
        #            f(x_res[i], ts[i])*h +
        #            f2(f(x_res[i], ts[i]), ts[i]) * 1.0/2.0 * h^2
        x_res[i+1] = x_res[i] +
                    f(x_res[i], ts[i])*h +
                    f2(f(x_res[i], ts[i]), ts[i]) * (h^2.0)/2.0
    end

    return x_res
end

function expl_euler_vec_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros((size(ts,1)+1, size(x_start,1)))
    x_res[1,:] = x_start

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res[i+1,:] = x_res[i,:] + f(x_res[i,:], ts[i])*h
    end

    return x_res
end

"""
Explicit Runge Kutta with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance
l      ... number of steps
bs     ... bs from the tableau
cs     ... cs from the tableau
as     ... as from the tableau

return a vector with inbetween steps
"""
function expl_runge_kutta_all(t_start, t_end, x_start, f, h, l, bs, cs, as)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros((size(ts,1)+1, size(x_start,1)))
    x_res[1] = x_start

    for i = 1:(size(ts,1))
        ks = zeros(l)
        ks[1] = f(x_res[i], ts[i])
        for j in 2:l
            ks[j] = f(x_res[i] .+ h*sum(as[j, 1:(j-1)].*ks[1:(j-1)]) ,ts[i] + cs[j]*h)
        end

        x_res[i+1] = x_res[i] .+ h*sum(bs.*ks)
    end

    return x_res
end

"""
Standard fourth order Runge Kutta with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return a vector with inbetween steps
"""
function standard_runge_kutta_all(t_start, t_end, x_start, f, h)
    bs = [1/6, 1/3, 1/3, 1/6]
    cs = [0, 1/2, 1/2, 1]
    as = [0.0  0.0  0.0  0.0;
          0.5  0.0  0.0  0.0;
          0.0  0.5  0.0  0.0;
          0.0  0.0  1.0  0.0]

    return expl_runge_kutta_all(t_start, t_end, x_start, f, h, 4, bs, cs, as)
end

"""
Modified euler (second order runge kutta)
"""
function second_order_runge_kutta_all(t_start, t_end, x_start, f, h)
    bs = [0, 1/2]
    cs = [0, 1]
    as = [0.0  0.0;
          0.5  0.0]

    return expl_runge_kutta_all(t_start, t_end, x_start, f, h, 2, bs, cs, as)
end

"""
Heun (second order runge kutta)
"""
function heun_runge_kutta_all(t_start, t_end, x_start, f, h)
    bs = [0, 1]
    cs = [0.5, 0.5]
    as = [0.0  0.0;
          1.0  0.0]

    return expl_runge_kutta_all(t_start, t_end, x_start, f, h, 2, bs, cs, as)
end

"""
Explicit Runge Kutta with inbetween steps for vectors
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance
l      ... number of steps
bs     ... bs from the tableau
cs     ... cs from the tableau
as     ... as from the tableau

return a vector with inbetween steps
"""
function expl_runge_kutta_vec_all(t_start, t_end, x_start, f, h, l, bs, cs, as)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros((size(ts,1)+1, size(x_start,1)))
    x_res[1,:] = x_start

    for i = 1:(size(ts,1))
        ks = zeros((l,size(x_start,1)))
        ks[1,:] = f(x_res[i,:], ts[i])

        for j in 2:l
            ks[j,:] = f(x_res[i,:] + h*sum(as[j, 1:(j-1)].*ks[1:(j-1), :], dims=1)' ,ts[i] + cs[j]*h)
        end

        x_res[i+1, :] = x_res[i,:] + h*sum(bs.*ks, dims=1)'
    end

    return x_res
end

"""
step for adaptive runge kutta
"""
function rk_step(x_j, t_j, h, l, as, bs, cs)

    ks = zeros((l,size(x_j,1)))
    ks[1,:] = f(x_j, t_j)

    for j in 2:l
        ks[j,:] = f(x_j + h*sum(as[j, 1:(j-1)].*ks[1:(j-1), :], dims=1)' ,t_j + cs[j]*h)
    end

    return x_j + h*sum(bs[1,:].*ks, dims=1)', x_j + h*sum(bs[2,:].*ks, dims=1)'

end

"""
Adaptive explicit Runge Kutta with inbetween steps for vectors
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h_start... initial sample point distance
l      ... number of steps
bs     ... bs from the tableau, matrix, first row is the higher order method
cs     ... cs from the tableau
as     ... as from the tableau
tol    ... tolerance of the adaptive method
min_order ... order of the worse method
h_min  ... minimal h where the algorithm will continue although not meeting the tolerance value

return a vector with inbetween steps and a vector with time steps
"""
function adaptive_expl_runge_kutta_vec_all(t_start, t_end, x_start, f, h_start, l, bs, cs, as, tol, min_order, break_on_min_h =true, h_min=0.00001)

    ts = [t_start]
    x_res = []
    push!(x_res, x_start)

    i = 1
    h = h_start

    while ts[end] < t_end
        #h = h_start

        # method with higher order, method with lower order
        x_1, xhat_1 = rk_step(x_res[end], ts[end], h, l, as, bs, cs)

        err = error_est(x_res[end], x_1, xhat_1)
        h = adaptive_h(tol, err, min_order, h)

        while err > tol && (h > h_min || !break_on_min_h)
            h_new = adaptive_h(tol, err, min_order, h)

            # no change, we got the minimal h
            if h_new == h
                break
            else
                h = h_new
            end

            x_1, xhat_1 = rk_step(x_res[end], ts[end], h, l, as, bs, cs)

            err = error_est(x_res[end], x_1, xhat_1)
        end

        if ts[i]+h > t_end
            h = t_end - ts[i]
        end

        push!(ts, ts[end]+h)
        push!(x_res, x_1)

        i += 1
    end

    return x_res, ts
end

"""
Standard fourth order Runge Kutta with inbetween steps for vectors
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return a vector with inbetween steps
"""
function standard_runge_kutta_vec_all(t_start, t_end, x_start, f, h)
    bs = [1/6, 1/3, 1/3, 1/6]
    cs = [0, 1/2, 1/2, 1]
    as = [0.0  0.0  0.0  0.0;
          0.5  0.0  0.0  0.0;
          0.0  0.5  0.0  0.0;
          0.0  0.0  1.0  0.0]

    return expl_runge_kutta_vec_all(t_start, t_end, x_start, f, h, 4, bs, cs, as)
end

"""
5-step Adam-Bashforth with inbetween steps for vector
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return a vector with inbetween steps
"""
function adams_bashforth_5_vec_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros((size(ts,1)+1, size(x_start,1)))
    x_res[1:5,:] = expl_euler_vec_all(t_start, t_start+h*4, x_start, f, h)

    for i = 5:(size(ts,1))
        x_res[i+1,:] = x_res[i,:] + h*(
                        (1901/720) * f(x_res[i,:], ts[i]) +
                        (-2774/720)* f(x_res[i-1,:], ts[i-1]) +
                        (2616/720) * f(x_res[i-2,:], ts[i-2]) +
                        (-1274/720)* f(x_res[i-3,:], ts[i-3]) +
                        (251/720)  * f(x_res[i-4,:], ts[i-4]))
    end

    return x_res
end

"""
k-step Adams Bashforth
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance
k      ... number of steps k in [1,6]
init_val_meth  ... one step method to obtain the first k steps
"""
function adams_bashforth_all(t_start, t_end, x_start, f, h, k, init_val_meth)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros(size(ts,1)+1)
    x_res[1:k] = init_val_meth(t_start, t_start+h*(k-1), x_start, f, h)

    coeffs = zeros(k)
    if k == 1
        coeffs[1] = 1
    elseif k == 2
        coeffs = [3, -1] * 1/2
    elseif k == 3
        coeffs = [23, -16, 5] * 1/12
    elseif k == 4
        coeffs = [55, -59, 37, 9] * 1/24
    elseif k == 5
        coeffs = [1901, -2774, 2616, -1274, 251] * 1/720
    elseif k == 6
        coeffs = [4277, -7923, 9982, -7298, 2877, -475] * 1/1440
    else
        return
    end

    for i = k:(size(ts,1))
        x_res[i+1] = x_res[i]
        # x_j+1 = x_j + h*coeff_0*x_j + h*coeff_1*x_j-1+....
        for j = 1:k
            x_res[i+1] += h*coeffs[j]* f(x_res[i-j+1], ts[i-j+1])
        end
    end

    return x_res
end

"""
k-step Nyström
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance
k      ... number of steps k in [1,6]
init_val_meth  ... one step method to obtain the first k steps
"""
function nystroem_all(t_start, t_end, x_start, f, h, k, init_val_meth)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros(size(ts,1)+1)
    x_res[1:k] = init_val_meth(t_start, t_start+h*(k-1), x_start, f, h)

    coeffs = zeros(k)
    if k == 2
        coeffs = [2, 0]
    elseif k == 3
        coeffs = [7, -2, 1] * 1/3
    elseif k == 4
        coeffs = [8, -5, 4, -1] * 1/3
    elseif k == 5
        coeffs = [269, -266, 294, -146, 29] * 1/90
    elseif k == 6
        coeffs = [297, -406, 574, -426, 169, -28] * 1/90
    else
        return
    end

    # x_j+1 = x_j-1 + h*coeff_0*x_j + h*coeff_1*x_j-1+....
    for i = k:(size(ts,1))
        x_res[i+1] = x_res[i-1]
        for j = 1:k
            x_res[i+1] += h*coeffs[j]* f(x_res[i-j+1], ts[i-j+1])
        end
    end

    return x_res
end

"""
k-step Adams Moulton
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... differential of f by x
h      ... sample point distance
k      ... number of steps k in [1,6]
init_val_meth  ... one step method to obtain the first k steps
"""
function adams_moulton_all(t_start, t_end, x_start, f, df, h, k, init_val_meth)

    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1:k] = init_val_meth(t_start, t_start+h*(k-1), x_start, f, h)

    coeffs = zeros(k)
    if k == 1
        coeffs[1] = 1
    elseif k == 2
        coeffs = [1, 1] * 1/2
    elseif k == 3
        coeffs = [5, 8, -1] * 1/12
    elseif k == 4
        coeffs = [9, 19, -5, 1] * 1/24
    elseif k == 5
        coeffs = [251, 646, -264, 106, -19] * 1/720
    elseif k == 6
        coeffs = [475, 1427, -798, 482, -173, 27] * 1/1440
    else
        return
    end

    for i = k:(size(ts,1)-1)
        # x_j+1 = x_j + h*coeff_0*x_j+1 + h*coeff_1*x_j+h*coeff_2*x_j-1+...
        x_approx = x_res[i]
        sum_bf = 0
        for l in 2:k
            sum_bf += h*coeffs[l]*f(x_res[i-l+2], ts[i-l+2])
        end

        for j = 1:10
            #use newton's method to approximate x_j+1
            # g(x_j+1) 0 = x_j+1 - x_j - h*b*f(x_j+1, t_j+1) + h*sum(b*f)
            # g'         = 1 - h*b*df
            #newton: a_n+1 = a_n - g(a_n)/g'(a_n)
            x_approx = x_approx -
                      (x_approx-x_res[i]-h*coeffs[1]*f(x_approx,ts[i+1])-sum_bf)/
                        (1-h*coeffs[1]*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_approx
    end

    return x_res
end

"""
k-step Milne Simpson
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... differential of f by x
h      ... sample point distance
k      ... number of steps k in [1,6]
init_val_meth  ... one step method to obtain the first k steps
"""
function milne_simpson_all(t_start, t_end, x_start, f, df, h, k, init_val_meth)

    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1:k] = init_val_meth(t_start, t_start+h*(k-1), x_start, f, h)

    coeffs = zeros(k)
    if k == 2
        coeffs = [1, 4, 1] * 1/3
    elseif k == 4
        coeffs = [29, 124, 24, 4, -1] * 1/90
    elseif k == 5
        coeffs = [28, 129, -14, 14, -6, 1] * 1/90
    else
        return
    end

    for i = k:(size(ts,1)-1)
        # x_j+1 = x_j-1 + h*coeff_0*x_j+1 + h*coeff_1*x_j+h*coeff_2*x_j-1+...
        x_approx = x_res[i]
        sum_bf = 0
        for l in 2:(k+1)
            sum_bf += h*coeffs[l]*f(x_res[i-l+2], ts[i-l+2])
        end

        for j = 1:10
            #use newton's method to approximate x_j+1
            # g(x_j+1) 0 = x_j+1 - x_j-1 - h*b*f(x_j+1, t_j+1) + h*sum(b*f)
            # g'         = 1 - h*b*df
            #newton: a_n+1 = a_n - g(a_n)/g'(a_n)
            x_approx = x_approx -
                      (x_approx-x_res[i-1]-h*coeffs[1]*f(x_approx,ts[i+1])-sum_bf)/
                        (1-h*coeffs[1]*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_res[i-1] + h*coeffs[1]*f(x_approx,ts[i+1]) + sum_bf
    end

    return x_res
end

"""
k-step BDF method
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... differential of f by x
h      ... sample point distance
k      ... number of steps k in [1,6]
init_val_meth  ... one step method to obtain the first k steps
"""
function BDF_all(t_start, t_end, x_start, f, df, h, k, init_val_meth)

    # coeff_0*x_j+1 + coeff_1*x_j + coeff_2 * x_j-1 +.... = h * f(x_j+1,t_j+1)

    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1:k] = init_val_meth(t_start, t_start+h*(k-1), x_start, f, h)

    coeffs = zeros(k)
    if k == 1
        coeffs = [1, -1]
    elseif k == 2
        coeffs = [1.5, -2.0, 0.5]
    elseif k == 3
        coeffs = [11.0/6, -3, 1.5, -1.0/3]
    elseif k == 4
        coeffs = [25.0/12, -4, 3, -4.0/3, 0.25]
    elseif k == 5
        coeffs = [137.0/60, -5, 5, -10.0/3, 1.25, -0.2]
    elseif k == 6
        coeffs = [147.0/60, -6, 15.0/2, -20.0/3, 15.0/4, -6.0/5, 1.0/6]
    else
        return
    end

    for i = k:(size(ts,1)-1)
        # coeff_1*x_j + coeff_2 * x_j-1 +....+coeffs_k*x_j+1-k
        x_approx = x_res[i]
        sum_prev = 0
        for l in 2:(k+1)
            sum_prev += coeffs[l]*f(x_res[i-l+2], ts[i-l+2])
        end

        for j = 1:10
            #use newton's method to approximate x_j+1
            # g(x_j+1) 0 = coeff_0*x_j+1 + coeff_1*x_j + ... - h*f(x_j+1, t_j+1)
            # g'         = 1 - h*df
            #newton: a_n+1 = a_n - g(a_n)/g'(a_n)
            x_approx = x_approx -
                      (coeffs[0]*x_approx + sum_prev - h*f(x_approx,ts[i+1]))/
                        (1-h*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_approx
    end

    return x_res
end

"""
Implicit Euler with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... differential of f by x
h      ... sample point distance
return a vector with inbetween steps
"""
function impl_euler_all(t_start, t_end, x_start, f, df, h, newton_steps=10)

    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1] = x_start

    for i = 1:(size(ts,1)-1)
        #implicit Euler
        #x_j+1 = x_j + f(t_j+1, x_j+1)*(t_j+1 - t_j)
        x_approx = x_res[i]
        for j = 1:newton_steps
            #use newton's method to approximate 0 = x_j+1 - x_j - h*f(x_j+1, t_j+1)
            #newton: a_n+1 = a_n - f(a_n)/f'(a_n)
            x_approx = x_approx - (x_approx-x_res[i]-h*f(x_approx,ts[i+1]))/(1-h*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_approx
    end

    return x_res
end

"""
Implicit Euler with inbetween steps for vectors
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... jacobian of f by x
h      ... sample point distance
return a vector with inbetween steps
"""
function impl_euler_vec_all(t_start, t_end, x_start, f, df, h, newton_steps)

    ts = collect(t_start:h:t_end)
    x_res = zeros((size(ts,1), size(x_start,1)))
    x_res[1,:] = x_start

    for i = 1:(size(ts,1)-1)
        neg_g(X) = -X + x_res[i,:] + h*f(X, ts[i+1])
        Dg(X) = LinearAlgebra.I(size(x_start,1)) - h*df(X, ts[i+1])

        x_res[i+1,:] = newton_method_vec(x_res[i,:], Dg, neg_g, newton_steps)
    end

    return x_res
end
