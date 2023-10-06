using StaticArrays
using LinearAlgebra
using DataStructures
using Statistics
using PyPlot

include("pattern.jl")

p3 = genPattern(π / 6, 3, SA[0.0, 0.0], r -> π / 6 .* @SVector ones(3))
p4 = genPattern(π / 4, 2, SA[0.0, 0.0])
p6 = genPattern(0.0, 3, SA[0.5, 0.0], r -> π .* @SVector ones(3))
pat = genPattern(0.0, 3, SVector{2}([0.0, 0.0]))

fermi(x) = 1 / (1 + exp(x))

function M(x, y, pat)
    r = SA[x, y]
    qr = pat.qs * r - pat.phases - pat.ψ(r)
    amp = 0.0
    if typeof(pat) <: HexPattern
        amp = cos(3.0 * mean(pat.ψ(r))) / 2.0
    end
    for i in 1:N(pat)
        @inbounds amp += cos(qr[i])
    end
    amp
end

function M(x, y, ps::Vector{<:AnalyticPattern})
    weight = fermi(y + 55)
    return sign(M(x, y, ps[1]) * weight + M(x, y, ps[2]) * (1.0 - weight))
end

function M3(x, y, x_cuts, ps)
    w0s = fermi.(x .- x_cuts)
    w1 = w0s[1]
    w2 = -w0s[1] + w0s[2]
    w3 = 1.0 - w0s[2]
    w1 * M(x, y, ps[1]) + w2 * M(x, y, ps[2]) + w3 * M(x, y, ps[3])
end

M0(x, y) = M(x, y, [p4, p6])
MC4(x, y) = M0(x, y)
M3(x, y) = M3(x, y, [-5, 5], [p3, p4, p6])
MC6(x, y) = M(x, y, p6)
MS6(x, y) = M(x, y, [p3, p6])

xs = -50:0.03:50
ys = -150.02:0.03:-50.03
