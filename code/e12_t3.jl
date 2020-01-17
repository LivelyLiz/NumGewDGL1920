euler = MathConstants.e
pi = MathConstants.pi

function alphas(k)
    if k == 1
        return [1, -1]
    elseif k == 2
        return [1.5, -2.0, 0.5]
    elseif k == 3
        return [11.0/6, -3, 1.5, -1.0/3]
    elseif k == 4
        return [25.0/12, -4, 3, -4.0/3, 0.25]
    elseif k == 5
        return [137.0/60, -5, 5, -10.0/3, 1.25, -0.2]
    elseif k == 6
        return [147.0/60, -6, 15.0/2, -20.0/3, 15.0/4, -6.0/5, 1.0/6]
    else
        return
    end
end

function root_locus(phi, coeffs)

    roh = 0
    sigma = (euler.^(1im * phi)).^(length(coeffs)-1)
    for i in 1:length(coeffs)
        roh = roh .+ coeffs[length(coeffs)-i+1]*(euler.^(1im*phi)).^(length(coeffs)-i+1)
    end

    return roh./sigma
end

phis = collect(0:0.01:2*pi)
plots = []
colors = ["red", "blue", "green", "magenta", "orange", "brown"]

for j in 1:6

    coeffs = alphas(j)
    res = root_locus(phis, coeffs)
    real_rl = real.(res)
    imag_rl = imag.(res)

    index = trunc(Int, length(real_rl)/2)

    subplot = plot(real_rl[1:index], imag_rl[1:index], xlabel="Re", ylabel = "Im", color=colors[j], arrow=true, aspect_ratio=1, title="k = $j")
    plot!(real_rl[index:(index*2)], imag_rl[index:(index*2)], arrow=true, color=colors[j])
    push!(plots, subplot)
end

p = plot(plots..., size=(1024,1024))

savefig(p, "../plot/e12_t3.pdf")
