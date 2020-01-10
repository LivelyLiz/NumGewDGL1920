using Plots
pyplot()

#Task 1

# a)
function R_a(z, s)
    res = 0
    for k = 0:s
        res += z^k/factorial(k)
    end
    return res
end

subplots_a = []

for i in 1:5
    subplot = contour(collect(-5:0.1:5), collect(-5:0.1:5), (x, y) -> abs(R_a(x + y*1im, i)), fill=true, levels=[0, 1], title = "ERK Order $i", xlabel = "Re", ylabel = "Im", aspect_ratio=1, legend=false)
    push!(subplots_a, subplot)
end

plot_a = plot(subplots_a...)

savefig(plot_a, "../plot/e11_t1_a.pdf")

# b)
R_b_1(z) = 1/12.0*z^4.0 + 1/6*z^3.0 + 1/2*z^2.0 + z + 1
R_b_2(z) = 3/18*z^3.0 + 1/2*z^2.0 + z + 1

plot_b = plot(contour(-5:0.1:5, -5:0.1:5, (x, y) -> abs(R_b_1(x + y*1im)), fill=true, levels=[0, 1], title = "1", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false),
              contour(-5:0.1:5, -5:0.1:5, (x, y) -> abs(R_b_2(x + y*1im)), fill=true, levels=[0, 1], title = "2", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false))

savefig(plot_b, "../plot/e11_t1_b.pdf")

# d)
R_c_implRK(z) = 1/(1-z)
R_c_implMid(z) = (z+2)/(2-z)
R_c_SDIRK(z, gamma) = 1 + z/((1-gamma*z)^2.0) * (1-2*gamma*z+1/2*z)

gamma1 = (3 + sqrt(3))/6
gamma2 = (3 - sqrt(3))/6

plot_c = plot(contour(-5:0.1:5, -5:0.1:5, (x, y) -> abs(R_c_implRK(x + y*1im)), fill=true, levels=[0, 1], title = "Impl. Euler", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false),
              contour(-5:0.1:5, -5:0.1:5, (x, y) -> abs(R_c_implMid(x + y*1im)), fill=true, levels=[0, 1], title = "Impl. Midpoint", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false),
              contour(-1:0.1:12, -7:0.1:7, (x, y) -> abs(R_c_SDIRK(x + y*1im, gamma1)), fill=true, levels=[0, 1], title = "SDIRK 3, +", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false),
              contour(-15:0.1:1, -8:0.1:8, (x, y) -> abs(R_c_SDIRK(x + y*1im, gamma2)), fill=true, levels=[0, 1], title = "SDIRK 3, -", xlabel = "Re", ylabel = "Im", aspect_ratio = 1, legend=false),
              size=(1024, 1024))

savefig(plot_c, "../plot/e11_t1_c.pdf")
