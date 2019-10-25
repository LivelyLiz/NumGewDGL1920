"""
Explicit Euler
t_start... given start t value
t_end  ... point where we want our x(t) value
x_start... x(t_start)
f      ... f(t, x) differential function
h      ... sample point distance
"""
function expl_euler(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:t_end)
    x_res = x_start

    for i = 1:size(ts,1)
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res = x_res + f(x_res)*h
    end

    return x_res
end
