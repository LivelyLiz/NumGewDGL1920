#Task 4
include("./integration_meth.jl")
using Plots

#initialization
x_start = 0
x_end = 1
y_start = 1
h = 0.2

#functions
f(y, x) = y^(3/2)
y(x, c) = 4 * 1 ./ (x .+ c) .^ 2

xs = collect(x_start:h:x_end)
yhs = expl_euler_all(x_start, x_end, y_start, f, h)
cs = []

#calculate analytical solution cs for each yh
for (i,yh) in enumerate(yhs)
    c1 = -xs[i] + 2/sqrt(yh)
    c2 = -xs[i] - 2/sqrt(yh)
    push!(cs, [c1, c2])
end

#plot
plotly()
p_c1 = plot(xs, yhs, label="explicit euler", title="positive c", xticks=x_start:h:x_end, linewidth=3)
for i in 1:5
    plot!(p_c1, xs[i:end], y(xs[i:end], cs[i][1]), label="i="*string(i-1)*" c="*string(cs[i][1]))
end

print(yhs)
print("\n")

p_c2 = plot(xs, yhs, label="explicit euler", title="negative c", xticks=x_start:h:x_end, linewidth=3)
for i in 1:5
    ys = y(xs[i:end], cs[i][2])

    print("i = "*string(i-1)*"\n")
    print(ys)
    print("\n")
    plot!(p_c2, xs[i:end], ys, label="i="*string(i-1)*" c="*string(cs[i][2]))
end

p = plot(p_c1, p_c2, size=(1500, 500))
savefig(p, "../plot/e3_t4.html")
display(p)

pyplot()
p_c1 = plot(xs, yhs, label="explicit euler", title="positive c", xticks=x_start:h:x_end, linewidth=3)
for i in 1:5
    plot!(p_c1, xs, y(xs, cs[i][1]), label="i="*string(i-1)*" c="*string(cs[i][1]))
end

p_c2 = plot(xs, yhs, label="explicit euler", title="negative c", xticks=x_start:h:x_end, linewidth=3)
for i in 1:5
    plot!(p_c2, xs, y(xs, cs[i][2]), label="i="*string(i-1)*" c="*string(cs[i][2]))
end

p = plot(p_c1, p_c2, size=(1500, 500))
savefig(p, "../plot/e3_t4.pdf")
