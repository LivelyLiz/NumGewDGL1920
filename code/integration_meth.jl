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
Implicit Euler with inbetween steps
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
df     ... differential of f by y
h      ... sample point distance
return a vector with inbetween steps
"""
function impl_euler_all(t_start, t_end, x_start, f, df, h)

    ts = collect(t_start:h:t_end)
    x_res = zeros(size(ts,1))
    x_res[1] = x_start

    for i = 1:(size(ts,1)-1)
        #implicit Euler
        #x_j+1 = x_j + f(t_j+1, x_j+1)*(t_j+1 - t_j)
        x_approx = x_res[i]
        for j = 1:10
            #use newton's method to approximate 0 = x_j+1 - x_j - h*f(x_j+1, t_j+1)
            #newton: a_n+1 = a_n - f(a_n)/f'(a_n)
            x_approx = x_approx - (x_approx-x_res[i]-h*f(x_approx,ts[i+1]))/(1-h*df(x_approx,ts[i+1]))
        end
        x_res[i+1] = x_approx
    end

    return x_res
end
