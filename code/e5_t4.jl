#Task 4
# see also e5_t2

function runge_kutta_vec_all(t_start, t_end, x_start, f, h)

    ts = collect(t_start:h:(t_end-h))
    x_res = zeros((size(ts,1)+1, size(x_start,1)))
    x_res[1,:] = x_start

    k1(xj, tj) = f(xj, tj)
    k2(xj, tj) = f(xj + h*1/2*k1(xj,tj), tj + 1/2*h)
    k3(xj, tj) = f(xj + h*k2(xj, tj), tj + h)
    k4(xj, tj) = f(xj + h*k3(xj, tj), tj + h)

    for i = 1:(size(ts,1))
        #explicit Euler
        #x_j+1 = x_j + f(t_j, x_j)*(t_j+1 - t_j)
        x_res[i+1,:] = x_res[i,:] + h*(1/6*k1(x_res[i,:], ts[i]) + 2/3 * k2(x_res[i,:],ts[i]) + 1/6*k4(x_res[i,:],ts[i]))
    end

    return x_res
end
