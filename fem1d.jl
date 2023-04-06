using Interpolations
# using QuadGK
using Plots

const Grid = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

function fem1d_lin_simple(f::Function, grid::Grid)
    M = size(grid)[1] - 2
    println(M)
    A = zeros(M, M)
    # b = (i -> quadgk(x -> phi(i-1, grid)(x) * f(x), 0, 1, rtol=1e-8)[1]).(1:M)
    b = (i -> num_int(x -> phi(i, grid)(x) * f(x), 0, 1, 1000)[1]).(1:M) # (phi, f)
    # println(b)

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

    # function a(u::Function, v::Function, x_l::Float64, x_r::Float64)
    #     num_int_simpsons(x -> p(x)*num_diff(u, x, x_l, x_r)*num_diff(v, x, x_l, x_r) + q(x)*u(x)*v(x), x_l, x_r, perc)
    #     # num_int(x -> p(x)*num_diff(u, x, x_l, x_r)*num_diff(v, x, x_l, x_r) + q(x)*u(x)*v(x), x_l, x_r, perc)
    # end

    x_l = grid[1]
    # println(x_l)
    x_r = grid[size(grid)[1]]
    # println(x_r)

    M = size(grid)[1] - 1
    # println(M)
    A = zeros(M+1, M+1)
    # b = (i -> quadgk(x -> phi(i-1, grid)(x) * f(x), 0, 1, rtol=1e-8)[1]).(1:M)
    # b = (i -> num_int(x -> phi(i, grid)(x) * f(x), x_l, x_r, 1000)[1]).(1:M) # (phi, f)
    # b = (i -> num_int_simpsons(x -> phi(i, grid)(x) * f(x), x_l, x_r, perc)[1]).(1:M) # (phi, f)
    b = zeros(M+1)
    # println(b)

    # h = (i -> grid[i+1]-grid[i]).(1:M+1)

    # for i=1:M
    #     tmp = num_int_simpsons(x -> p(x) * phi_prime_right(i, grid)^2 + q(x) * (phi(i, grid)(x))^2, grid[i], grid[i+1], perc)
    #     # tmp2 = num_int_simpsons(x -> p(x) * phi_prime_left(i+1, grid)^2 + q(x) * (phi(i+1, grid)(x))^2, grid[i], grid[i+1], perc)
    #     tmp2 = num_int_simpsons(x -> p(x) * phi_prime_right(i, grid) * phi_prime_left(i+1, grid) + q(x) * phi(i, grid)(x) * phi(i+1, grid)(x), grid[i], grid[i+1], perc)
    #     # println(tmp)
    #     # println(tmp2)
    # end
    A[1, 1] += num_int_simpsons(x -> p(x) * phi_prime_right(1, grid)^2 + q(x) * (phi(1, grid)(x))^2, grid[1], grid[1+1], perc)
    b[1] += num_int_simpsons(x -> phi(1, grid)(x) * f(x), grid[1], grid[1+1], perc)
    # tmp = num_int_simpsons(x -> (phi(1, grid)(x))^2, -0.2, -0.1, perc)
    for i=1:M
        # A[i, i] += a(x -> phi(i, grid)(x), x -> phi(i, grid)(x), x_l, x_r)
        # A[i, i+1] = a(x -> phi(i, grid)(x), x -> phi(i+1, grid)(x), x_l, x_r)
        # A[i+1, i] = A[i, i+1]
        # A[i+1, i+1] += a(x -> phi(i+1, grid)(x), x -> phi(i+1, grid)(x), x_l, x_r)

        A[i, i] += num_int_simpsons(x -> p(x) * phi_prime_right(i, grid)^2 + q(x) * (phi(i, grid)(x))^2, grid[i], grid[i+1], perc)
        A[i+1, i] = num_int_simpsons(x -> p(x) * phi_prime_right(i, grid) * phi_prime_left(i+1, grid) + q(x) * phi(i, grid)(x) * phi(i+1, grid)(x), grid[i], grid[i+1], perc)
        A[i, i+1] = A[i+1, i]
        A[i+1, i+1] += num_int_simpsons(x -> p(x) * phi_prime_left(i+1, grid)^2 + q(x) * (phi(i+1, grid)(x))^2, grid[i], grid[i+1], perc)

        b[i] += num_int_simpsons(x -> phi(i, grid)(x) * f(x), grid[i], grid[i+1], perc)
        b[i+1] += num_int_simpsons(x -> phi(i+1, grid)(x) * f(x), grid[i], grid[i+1], perc)
    end
    b[M+1] += num_int_simpsons(x -> phi(M, grid)(x) * f(x), grid[M], grid[M+1], perc)

    A[M+1, M+1] += num_int_simpsons(x -> p(x) * (phi_prime_right(M, grid))^2 + q(x) * phi(M, grid)(x)^2, grid[M], grid[M+1], perc)

    # A[M, M] = num_int(x -> p(x) * (phi_prime(M, grid, x))^2 + q(x) * phi(M, grid)(x)^2, x_l, x_r, 1000)
    # A[M, M] = a(x -> phi(M, grid)(x), x -> phi(M, grid)(x), x_l, x_r)
    # A[M, M] = 30
    # println(A[M, M])

    # println(phi_prime(1, grid, 0.05))
    # println(grid[1])

    # A = A[2:M, 2:M]
    # b = b[2:M]
    # A = A[1:M, 1:M]
    # A[1,1] = A[1,1] * 200
    # b = b[1:M]

    A[1,1] = 1
    A[1,2] = 0
    # A[2,1] = 0
    b[1] = 0

    A = A[1:M, 1:M]
    b = b[1:M]

    display(b)
    display(A)

    res = A\b
    # res[1] = 0
    # res[M+1] = 0
    # prepend!(res, 0)
    # append!(res, 0)
    display(res)
    res
end


function phi(i::Integer, grid::Grid)
    xs = (x -> round((x-grid[1])*(1/(grid[2]-grid[1]))) == (i-1) ? 1 : 0).(grid)
    linear_interpolation(grid, xs)
end

function phi_prime_left(i, grid)
    return 1 / (grid[i] - grid[i-1])
    # return phi_prime((grid[i] + grid[i-1])/2, i, grid)
end

function phi_prime_right(i, grid)
    return -1 / (grid[i+1] - grid[i])
    # return phi_prime((grid[i] + grid[i+1])/2, i, grid)
end

function phi_prime(x, i, grid)
    return num_diff(k -> phi(i, grid)(k), x, -1.0, 1.0)
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
    # s += 4 * sum(f((1:2:n) .* h .+ low))
    for i=1:2:n
        s += 4 * f(i * h + low)
    end
    for i=2:2:n-1
        s += 2 * f(i * h + low)
    end
    # s += 2 * sum(f((2:2:n-1) .* h .+ low))
    return (h/3) * s
    # for i = low:1/per:(high-(1/per))
    #     a = i
    #     b = max(i + 1/per, high)
    #     sum += ((b-a)/6) * (f(a) + 4*f((a+b)/2) + f(b))
    # end
    # sum
end

# u(x) = (x * (1-x)) / 2
# f(x) = 1

# u(x) = sin(pi * x)
# f(x) = pi^2 * sin(pi*x)

# begin
#     u(x) = sin(2 * pi * x)
#     f(x) = 4 * pi^2 * sin(2 * pi*x)

#     h = 100
#     grid = 0:1/h:1

#     xi = fem1d_lin_simple(f, grid)
#     plot(x -> u_h(x, xi, grid), 0:0.001:1)
#     plot!(u, 0:0.001:1)
# end

begin
    c = 5
    p(x) = x^2 + 1
    q(x) = c - x^2
    # f(x) = x^4 - (7 + c)*(x^2) + (c - 2)
    f(x) = 2 + c + 5*x^2 - c*x^2 + x^4
    u(x) = 1-x^2

    h1 = 1 / 100
    h2 = h1 / 2
    h3 = h1 / 10
    grid = -1:h1:1
    grid2 = -1:h2:1

    xi = fem1d_lin_sl(f, p, q, grid)
    xi2 = fem1d_lin_sl(f, p, q, grid2)
    plot(x -> u_h(x, xi, grid), -1:h3:1)
    plot!(x -> u_h(x, xi2, grid2), -1:h3:1)
    plot!(u, -1:h3:1)

    d_inf1 = 0
    d_inf2 = 0
    for i=-1:h3:1
        global d_inf1 = max(d_inf1, abs(u(i) - u_h(i, xi, grid)))
        global d_inf2 = max(d_inf2, abs(u(i) - u_h(i, xi2, grid2)))
    end

    r = log(d_inf2 / d_inf1) / log(h2 / h1)
    println("r is")
    println(d_inf1)
    println(d_inf2)
    println(r)
end
