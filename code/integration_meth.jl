"""
Explicit Euler
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(x, t) differential function
h      ... sample point distance

return the value for t_end
"""
function expl_euler(t_start, t_end, x_start, f, h)

    ts = collect(t_start+h:h:t_end)
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
h      ... sample point distance

return a vector with inbetween steps
"""
function expl_euler_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start+h:h:t_end)
    x_res = zeros(size(ts,1)+1)
    x_res[1] = x_start

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res[i+1] = x_res[i] + f(x_res[i], ts[i])*h
    end

    return x_res
end

function expl_euler_vec_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start+h:h:t_end)
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
