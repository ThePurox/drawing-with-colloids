using LinearAlgebra
using StaticArrays
using Statistics
using AdaptiveBD
using PyPlot
using Folds
using ProgressMeter

include("plot-tools.jl")
include("tools.jl")
include("forces.jl")
include("gen-patterns.jl")
include("loops.jl")

const north = 0.0
const south = Float64(π)

gen_mag(m_func) = m_func.(xx, yy)

gen_mags() = Folds.map(gen_mag, m_funcs)
ip(mag) = InterpolatedPattern(mag=mag, xs=xs, ys=ys, z=0.5)
gen_cpats(mags) = Folds.map(ip, mags)

function d(a, t)
    loop4a(a,
           SA[(1.0 - dϕi, north), (1.0 + dϕi, south), (2.0 - dϕi, north),
              (2.0 + dϕi, south)], t)
end # down
function u(a, t)
    loop4a(a,
           SA[(3.0 - dϕi, north), (3.0 + dϕi, south), (4.0 - dϕi, north),
              (4.0 + dϕi, south)], t)
end # up
function l(a, t)
    loop4a(a,
           SA[(0.0 - dϕi, north), (0.0 + dϕi, south), (1.0 - dϕi, north),
              (1.0 + dϕi, south)], t)
end # left
function r(a, t)
    loop4a(a,
           SA[(2.0 - dϕi, north), (2.0 + dϕi, south), (3.0 - dϕi, north),
              (3.0 + dϕi, south)], t)
end # right

lu(a, t) = loop4a(a, SA[(0.0 - dϕi, north), (0.0 + dϕi, south)], t) # left up
rd(a, t) = loop4a(a, SA[(2.0 - dϕi, north), (2.0 + dϕi, south)], t) # right down
ld(a, t) = loop4a(a, SA[(1.0 - dϕi, north), (1.0 + dϕi, south)], t) # left down
ru(a, t) = loop4a(a, SA[(3.0 - dϕi, north), (3.0 + dϕi, south)], t) # right up

function luu(a, t)
    loop4a(a,
           SA[(0.0 - dϕi, north), (0.0 + dϕi, south), (3.0 - dϕi, north),
              (3.0 + dϕi, south), (4.0 - dϕi, north), (4.0 + dϕi, south)], t)
end # left up
function rdd(a, t)
    loop4a(a,
           SA[(2.0 - dϕi, north), (2.0 + dϕi, south), (1.0 - dϕi, north),
              (1.0 + dϕi, south), (2.0 - dϕi, north), (2.0 + dϕi, south)], t)
end # right down
function ldd(a, t)
    loop4a(a,
           SA[(1.0 - dϕi, north), (1.0 + dϕi, south), (1.0 - dϕi, north),
              (1.0 + dϕi, south), (2.0 - dϕi, north), (2.0 + dϕi, south)], t)
end # left down
function ruu(a, t)
    loop4a(a,
           SA[(3.0 - dϕi, north), (3.0 + dϕi, south), (3.0 - dϕi, north),
              (3.0 + dϕi, south), (4.0 - dϕi, north), (4.0 + dϕi, south)], t)
end # right up

loop4a(α, loop, t) = loop4p(0.5 - α, loop, t)

include("alphabet-loops.jl")

function write_alphabet_shift(t)
    n_all = floor(Int64, t / (Nₚ * T()))
    t0 = T() * n_all
    ta = t % T() + t0
    n_inner = floor(Int, ((t / T()) % Nₚ))
    if n_inner == 0
        write_A((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 1
        write_B((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 2
        write_C((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 3
        write_D((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 4
        write_E((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 5
        write_F((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 6
        write_G((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 7
        write_H((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 8
        write_I((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 9
        write_J((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 10
        write_K((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 11
        write_L((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 12
        write_M((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 13
        write_N((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 14
        write_O((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 15
        write_P((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 16
        write_Q((2n_inner + 1) * dϕ, ta)
    elseif n_inner == 17
        write_R((2n_inner + 1) * dϕ, ta)
    else
        return neut
    end
end

function write_alphabet(t)
    n_all = floor(Int64, t / (Nₚ * T()))
    t0 = T() * n_all
    ta = t % T() + t0
    n_inner = floor(Int, ((t / T()) % Nₚ))
    if n_inner == 0
        write_A((2n_inner) * dϕ, ta)
    elseif n_inner == 1
        write_B((2n_inner) * dϕ, ta)
    elseif n_inner == 2
        write_C((2n_inner) * dϕ, ta)
    elseif n_inner == 3
        write_D((2n_inner) * dϕ, ta)
    elseif n_inner == 4
        write_E((2n_inner) * dϕ, ta)
    elseif n_inner == 5
        write_F((2n_inner) * dϕ, ta)
    elseif n_inner == 6
        write_G((2n_inner) * dϕ, ta)
    elseif n_inner == 7
        write_H((2n_inner) * dϕ, ta)
    elseif n_inner == 8
        write_I((2n_inner) * dϕ, ta)
    elseif n_inner == 9
        write_J((2n_inner) * dϕ, ta)
    elseif n_inner == 10
        write_K((2n_inner) * dϕ, ta)
    elseif n_inner == 11
        write_L((2n_inner) * dϕ, ta)
    elseif n_inner == 12
        write_M((2n_inner) * dϕ, ta)
    elseif n_inner == 13
        write_N((2n_inner) * dϕ, ta)
    elseif n_inner == 14
        write_O((2n_inner) * dϕ, ta)
    elseif n_inner == 15
        write_P((2n_inner) * dϕ, ta)
    elseif n_inner == 16
        write_Q((2n_inner) * dϕ, ta)
    elseif n_inner == 17
        write_R((2n_inner) * dϕ, ta)
    else
        return neut
    end
end

struct Loop
    T::Float64
    crossings::Vector{Float64}
end

function control_space(t, L::Loop)
    i = floor(Int, t / L.T) + 1
    t = (t % L.T) / L.T
    phi = L.crossings[i]
    if i % 2 == 1
        theta = π * t
    else
        theta = π * (1 - t)
    end
    phi, theta
end

cc = [write_alphabet.(T() .* t)[1]
      for t in 0:(1 / 16):2525 if abs(write_alphabet(T() * t)[2] - π / 2) < 1e-12]
lop = Loop(T(), cc)
write_alphabet_compact(t) = control_space(t, lop)

function plot_all_letters(datas)
    fig, axs = subplots(3, 6, sharex=true)
    for (i, d) in enumerate(datas)
        axs[i].set_aspect("equal")
        plot_pats(d[2], ax=axs[i])
    end
end
