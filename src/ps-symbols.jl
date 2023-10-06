using StaticArrays
using LinearAlgebra
using DataStructures
using Statistics
using AdaptiveBD
using Base.Threads
using PyPlot
using Folds
using ProgressMeter

include("plot-tools.jl")
include("tools.jl")
include("forces.jl")
include("gen-patterns.jl")
include("loops.jl")

const Lsq = 20
const Lx = 20
const LΔ = 26

function gen_mags()
    p1 = genPattern(π / 4 + 0.0, 2, SA[0.0, 0.0])
    p2 = genPattern(π / 4 + π / 6, 2, SA[0.0, 0.0])
    p3 = genPattern(π / 4 + π / 3, 2, SA[0.0, 0.0])

    p6 = genPattern(0.0, 3, SA[0.5, 0.0], r -> π * SA[1, 1, 1])

    ψs(r) = SA[0.0, 0.0, atan(r[1], r[2])]
    pat = genPattern(0.0, 3, SVector{2}([0.0, 0.0]), ψs)

    m1(x, y) = M(x, y, [p1, p6])
    m2(x, y) = M(x, y, [p2, p6])
    m3(x, y) = M(x, y, [p3, p6])

    xx, yy = meshgrid(xs, ys)

    m41 = sign.(m1.(xx, yy))
    m42 = sign.(m2.(xx, yy))
    m43 = sign.(m3.(xx, yy))
    mag50 = sign.([M(SA[y, x], pat) for x in xs, y in xs])

    mag1 = vcat(m41, mag50)
    mag2 = vcat(m42, mag50)
    mag3 = vcat(m43, mag50)
    return mag1, mag2, mag3
end

mm1(x, y) = M(x, y, p1)
mm2(x, y) = M(x, y, p2)
mm3(x, y) = M(x, y, p3)

loop41(loop, t) = loop4p(0.5, loop, t)
loop42(loop, t) = loop4p(0.5 + 1 / 3, loop, t)
loop43(loop, t) = loop4p(0.5 + 2 / 3, loop, t)

function write_square(t)
    n = 0.0
    s = π
    # dϕ = 0.05
    dϕ = 0.15
    if t < Lsq / 2 * T()
        loop41(SA[(-dϕ, n), (dϕ, s), (1.0 - dϕ, n), (1.0 + dϕ, s)], t) # left
    elseif t < 3Lsq / 2 * T()
        loop41(SA[(1.0 - dϕ, n), (1.0 + dϕ, s), (2.0 - dϕ, n), (2.0 + dϕ, s)], t) # down
    elseif t < 5Lsq / 2 * T()
        loop41(SA[(2.0 - dϕ, n), (2.0 + dϕ, s), (3.0 - dϕ, n), (3.0 + dϕ, s)], t) # right
    elseif t < 7Lsq / 2 * T()
        loop41(SA[(3.0 - dϕ, n), (3.0 + dϕ, s), (4.0 - dϕ, n), (4.0 + dϕ, s)], t) # up
    elseif t < 8Lsq / 2 * T()
        loop41(SA[(-dϕ, n), (dϕ, s), (1.0 - dϕ, n), (1.0 + dϕ, s)], t) # left
    else
        loop4(SA[(0.0, n), (0.0, n)], t)
    end
end

function write_x(t)
    n = 0.0
    s = Float64(π)
    dϕ = 0.15
    if t < Lx / 2 * T()
        loop42(SA[(-dϕ, n), (dϕ, s)], t) # left up
    elseif t < 3Lx / 2 * T()
        loop42(SA[(2.0 - dϕ, n), (2.0 + dϕ, s)], t) # right down
    elseif t < 4Lx / 2 * T()
        loop42(SA[(-dϕ, n), (dϕ, s)], t) # left up
    elseif t < 5Lx / 2 * T()
        loop42(SA[(1.0 - dϕ, n), (1.0 + dϕ, s)], t) # right down
    elseif t < 7Lx / 2 * T()
        loop42(SA[(3.0 - dϕ, n), (3.0 + dϕ, s)], t) # right down
    elseif t < 8Lx / 2 * T()
        loop42(SA[(1.0 - dϕ, n), (1.0 + dϕ, s)], t) # right down
    else
        loop4(SA[(0.0, n), (0.0, n)], t)
    end
end

function write_Δ(t)
    n = 0.0
    s = Float64(π)
    dϕ = 0.15
    if t < 2LΔ / 2 * T()
        loop43(SA[(1.0 - dϕ, n), (1.0 + dϕ, s)], t) # right down
    elseif t < 4LΔ / 2 * T()
        loop43(SA[(2.0 - dϕ, n), (2.0 + dϕ, s), (3.0 - dϕ, n), (3.0 + dϕ, s)], t) # right
    elseif t < 6LΔ / 2 * T()
        loop43(SA[(-dϕ, n), (dϕ, s)], t) # left up
    else
        loop4(SA[(0.0, n), (0.0, n)], t)
    end
end

function write_symbols(t)
    t0 = T() * floor(Int64, t / 3T())
    ta = t % 3T() + t0
    if (t % 3T()) / T() < 1
        return write_square(ta)
    elseif (t % 3T()) / T() < 2
        return write_x(ta)
    elseif (t % 3T()) / T() < 3
        return write_Δ(ta)
    end
end

function ps_loop(t)
    nᵢₙ = 200
    height = 100

    nsq = Lsq * 4
    nx = Lx * 4
    nΔ = LΔ * 3
    nsymbols = 10max(LΔ, Lsq, Lx)

    if t < nᵢₙ * T()
        return spiral_in(t)
    elseif t < (nᵢₙ + height / 4 + 8) * T()
        return loop6(SA[(3, 0), (10, Float64(π)), (0, 0), (8, Float64(π))], t)
    elseif t < (nᵢₙ + height / 2) * T()
        return loop6(SA[(4, 0), (10, Float64(π)), (0, 0), (8, Float64(π))], t)
    else
        write_symbols(3 * (t - (nᵢₙ + height / 2) * T()))
    end
end
