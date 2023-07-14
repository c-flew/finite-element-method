using Interpolations
using Plots

const Grid = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

function fem1d_lin_simple(f::Function, grid::Grid)
    M = size(grid)[1] - 2
    println(M)
    A = zeros(M, M)
    b = (i -> num_int(x -> phi(i, grid)(x) * f(x), 0, 1, 1000)[1]).(1:M) # (phi, f)

    h = (i -> grid[i+1]-grid[i]).(1:M+1)

    A[1, 1] = (1/h[1]) + (1/h[2])

    for i=2:M
        A[i, i] = (1/h[i]) + 1/(h[i+1]) # main diagonal
        A[i, i-1] = -1/h[i] # bottom adj diagonal
        A[i-1, i] = A[i, i-1] # top adj diagonal
    end

    display(A)

    A\b
end

function fem1d_lin_sl(f::Function, p::Function, q::Function, grid::Grid)
    perc = 1000

    x_l = grid[1]
    x_r = grid[size(grid)[1]]

    M = size(grid)[1] - 1
    A = zeros(M+1, M+1)
    b = zeros(M+1)

    A[1, 1] += quadrature(x -> p(x) * phi_prime_right(1, grid)^2 + q(x) * (phi(1, grid)(x))^2, grid[1], grid[1+1])
    b[1] += num_int_simpsons(x -> phi(1, grid)(x) * f(x), grid[1], grid[1+1], perc)

    for i=1:M
        A[i, i] += quadrature(x -> p(x) * phi_prime_right(i, grid)^2 + q(x) * (phi(i, grid)(x))^2, grid[i], grid[i+1])
        A[i+1, i] = quadrature(x -> p(x) * phi_prime_right(i, grid) * phi_prime_left(i+1, grid) + q(x) * phi(i, grid)(x) * phi(i+1, grid)(x), grid[i], grid[i+1])
        A[i, i+1] = A[i+1, i]
        A[i+1, i+1] += quadrature(x -> p(x) * phi_prime_left(i+1, grid)^2 + q(x) * (phi(i+1, grid)(x))^2, grid[i], grid[i+1])

        b[i] += num_int_simpsons(x -> phi(i, grid)(x) * f(x), grid[i], grid[i+1], perc)
        b[i+1] += num_int_simpsons(x -> phi(i+1, grid)(x) * f(x), grid[i], grid[i+1], perc)
    end
    b[M+1] += num_int_simpsons(x -> phi(M, grid)(x) * f(x), grid[M], grid[M+1], perc)

    A[M+1, M+1] += quadrature(x -> p(x) * (phi_prime_right(M, grid))^2 + q(x) * phi(M, grid)(x)^2, grid[M], grid[M+1])


    A[1,1] = 1
    A[1,2] = 0
    b[1] = 0

    A = A[1:M, 1:M]
    b = b[1:M]

    display(b)
    display(A)

    res = A\b
    display(res)
    res
end


function phi(i::Integer, grid::Grid)
    xs = (x -> round((x-grid[1])*(1/(grid[2]-grid[1]))) == (i-1) ? 1 : 0).(grid)
    linear_interpolation(grid, xs)
end

function phi_prime_left(i, grid)
    return 1 / (grid[i] - grid[i-1])
end

function phi_prime_right(i, grid)
    return -1 / (grid[i+1] - grid[i])
end

function phi_prime(x, i, grid)
    return num_diff(k -> phi(i, grid)(k), x, -1.0, 1.0)
end

function quadrature(f, a, b)
    w = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891]
    ab = [0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640]

    sum = 0
    for i=1:size(w)[1]
        sum += w[i] * f(((b-a)/2) * ab[i] + (a+b)/2)
    end
    sum * (b-a)/2
end

function u_h(x, xi, grid)
    M = size(xi)[1]
    sum = 0
    for i=1:M
        sum += xi[i] * phi(i, grid)(x)
    end
    sum
end

function num_int(f, low, high, per)
    @assert low <= high
    sum = 0
    for i=low:1/per:high
        sum += f(i) * 1/per
    end
    sum
end

function num_diff(f::Function, x::Float64, x_l::Float64, x_r::Float64, h::Float64=0.0001)
    up = min(x+h, x_r)
    low = max(x-h, x_l)
    (f(up) - f(low)) / (up - low)
end


function num_int_simpsons(f::Function, low, high, n=6)
    @assert low <= high
    s = 0.0
    h = (high - low) / n
    s += f(low) + f(high)
    for i=1:2:n
        s += 4 * f(i * h + low)
    end
    for i=2:2:n-1
        s += 2 * f(i * h + low)
    end
    return (h/3) * s
end

# u(x) = (x * (1-x)) / 2
# f(x) = 1

# u(x) = sin(pi * x)
# f(x) = pi^2 * sin(pi*x)

# begin
#     # u(x) = sin(2 * pi * x)
#     # f(x) = 4 * pi^2 * sin(2 * pi*x)
#     f(x) = x


#     h = 100
#     grid = 0:1/h:1

#     xi = fem1d_lin_simple(f, grid)
#     plot(x -> u_h(x, xi, grid), 0:0.001:1)
#     # plot!(u, 0:0.001:1)
# end

begin
    c = 5
    p(x) = x^2 + 1
    q(x) = c - x^2
    # f(x) = x^4 - (7 + c)*(x^2) + (c - 2)
    f(x) = 2 + c + 5*x^2 - c*x^2 + x^4
    u(x) = 1-x^2

    #p(x) = x^3+5
    #q(x) = cos(x)
    #f(x) = 5*sin(10*x)

    # p(x) = 1
    # q(x) = 0
    # f(x) = 5 * (exp(x^3) - 1)
    # f(x) = 10*sin(10 * pi * x)

    h1 = 1 / 5
    h2 = h1 / 2
    h3 = h1 / 10
    grid = -1:h1:1
    grid2 = -1:h2:1

    xi = fem1d_lin_sl(f, p, q, grid)
    xi2 = fem1d_lin_sl(f, p, q, grid2)

    d_inf1 = 0
    d_inf2 = 0
    for i=-1:h3:1
        global d_inf1 = max(d_inf1, abs(u(i) - u_h(i, xi, grid)))
        global d_inf2 = max(d_inf2, abs(u(i) - u_h(i, xi2, grid2)))
    end

    max_u = 0
    for i=-1:h3:1
        global max_u = max(max_u, abs(u(i)))
    end

    d_inf1 /- max_u
    d_inf2 /- max_u


    # d_l11 = num_int_simpsons(x -> abs(u(x) - u_h(x, xi, grid)), -1, 1, 1/h3)
    # d_l12 = num_int_simpsons(x -> abs(u(x) - u_h(x, xi2, grid2)), -1, 1, 1/h3)

    # d_l11 = quadrature(x -> abs(u(x) - u_h(x, xi, grid)), -1, 1)
    # d_l12 = quadrature(x -> abs(u(x) - u_h(x, xi2, grid2)), -1, 1)

    d_l11 = 0
    d_l12 = 0
    l1_u = 0
    for i=1:size(grid)[1]-1
        global d_l11 += quadrature(x -> abs(u(x) - u_h(x, xi, grid)), grid[i], grid[i+1])
        global l1_u += quadrature(x -> abs(u(x)), grid[i], grid[i+1])
    end
    for i=1:size(grid2)[1]-1
        global d_l12 += quadrature(x -> abs(u(x) - u_h(x, xi2, grid2)), grid2[i], grid2[i+1])
    end

    d_l11 /= l1_u
    d_l12 /= l1_u

    # d_l21 = num_int_simpsons(x -> abs(u(x) - u_h(x, xi, grid))^2, -1, 1, 1/h3) ^ 0.5
    # d_l22 = num_int_simpsons(x -> abs(u(x) - u_h(x, xi2, grid2))^2, -1, 1, 1/h3) ^ 0.5
    d_l21 = 0
    d_l22 = 0
    l2_u = 0
    for i=1:size(grid)[1]-1
        global d_l21 += quadrature(x -> abs(u(x) - u_h(x, xi, grid))^2, grid[i], grid[i+1])
        global l2_u += quadrature(x -> abs(u(x))^2, grid[i], grid[i+1])
    end
    d_l21 = sqrt(d_l21)
    for i=1:size(grid2)[1]-1
        global d_l22 += quadrature(x -> abs(u(x) - u_h(x, xi2, grid2))^2, grid2[i], grid2[i+1])
    end
    d_l22 = sqrt(d_l22)

    d_l21 /= l2_u
    d_l22 /= l2_u

    r = log(d_inf2 / d_inf1) / log(h2 / h1)
    println("r is")

    r2 = log(d_l12 / d_l11) / log(h2 / h1)
    r3 = log(d_l22 / d_l21) / log(h2 / h1)
    # println(d_inf1)
    # println(d_inf2)
    println(r)
    println(r2)
    println(r3)
    println(d_inf1)
    # println(d_inf2)
    println(d_l11)
    println(d_l21)

    # println(d_l21)
    # println(d_l21t)

    # println(d_l22)
    # println(d_l22t)

    # println(quadrature(x -> abs(u(x)-u_h(x, xi, grid))^2, -1, 1)^0.5)
    # println(quadrature(x -> abs(u(x)-u_h(x, xi2, grid2))^2, -1, 1)^0.5)
    # ttt = quadrature(x -> abs(u(x) - u_h(x, xi, grid)), -1, 1)
    # d_ttt = num_int_simpsons(x -> abs(u(x) - u_h(x, xi, grid)), -1, 1, 1/(100*h3))
    # println(d_ttt)
    # println(ttt)
    # plot(x -> u_h(x, xi, grid), -1:h3:1)
    plot(x -> u_h(x, xi2, grid2), -1:h3:1, label="u_h")
    plot!(u, -1:h3:1, label="u")
    savefig("plot.svg")
end
