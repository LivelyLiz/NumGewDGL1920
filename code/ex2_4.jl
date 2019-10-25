#Task 4
using Plots
plotly()

y(x) = 4/15*x.^(5/2) - 19/15*x .+ 1 #.+ and .^ for pointwise
ddy(x) = sqrt.(x)                   # sqrt. pointwise sqrt
h = 0.01                            # sample point distance
xs = collect(0:h:1) #collect -> make range in brackts to array

# y values
ys_real = y(xs) #values for y(x)   at sample points
ddys = ddy(xs)  #values for y''(x) at sample points
ys = ddys       #values for solving system of lin eq
ys[1] = ys_real[1]
ys[end] = ys_real[end]

# setting up the system of linear equations
A = zeros(size(xs,1), size(xs,1)) #zero matrix of size of number of sample points
#we know the first and last values
A[1,1] = h^2
A[end, end] = h^2

# setting up the coefficients
for i = 2:(size(xs,1)-1)
    A[i,i-1] = 1
    A[i, i] = -2
    A[i, i+1] = 1
end
A = A.*1/(h^2)

ys_num = A\ys #solving the system of linear equations

err = abs.(ys_num-ys_real) #getting the error

p = plot(xs, err, title="Error of numerical method", label="error", xlabel="x", ylabel="abs(y_num-y_real)")
savefig(p, "e2_t4.html")
display(p)
