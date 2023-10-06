import FFTW
import Interpolations

using StaticArrays
using Statistics

include("tools.jl")

α_per(q, p, n) = 2 * atan(q / p * sin(π / n) / (1.0 + q / p * cos(π / n)))
α_perd(q, p, n) = α_per(q, p, n) * 180 / π
αₘ(k, n) = α_per(1, n * k + 1, n)
αₘd(k, n) = αₘ(k, n) * 180 / π

abstract type AbstractPattern end
abstract type AnalyticPattern <: AbstractPattern end

struct InterpolatedPattern{T} <: AbstractPattern
    ∇H::SVector{2, Float64}
    Hₚs::SVector{3, T}
end

struct SquarePattern{T} <: AnalyticPattern
    α::Float64
    qs::SMatrix{2, 2, Float64, 4}
    phases::SVector{2, Float64}
    ψ::T
end

struct HexPattern{T} <: AnalyticPattern
    α::Float64
    qs::SMatrix{3, 2, Float64, 6}
    phases::SVector{3, Float64}
    ψ::T
end

N(p::SquarePattern{T}) where {T} = 2
N(p::HexPattern{T}) where {T} = 3
q(p::SquarePattern{T}) where {T} = 2π
q(p::HexPattern{T}) where {T} = 2π / sin(π / 3)

function q(i, N, φ₀)
    f = q₀ = 2π
    if N == 3
        q₀ *= 1 / sin(π / 3)
        f = 4π
    end
    φ = f * i / 2N + φ₀
    q₀ * SA[-sin(φ), cos(φ)]
end

q(N, φ₀) =
    if N == 2
        return [q(1, N, φ₀) q(2, N, φ₀)]
    elseif N == 3
        return [q(1, N, φ₀) q(2, N, φ₀) q(3, N, φ₀)]'
    end

function genPattern(α::Float64, N::Int64, shift::SVector{2}, ψ=nothing)
    qs = q(N, α)
    R = SA[cos(α) -sin(α); sin(α) cos(α)]
    shift = R * shift
    phases = qs * shift
    if N == 3
        if ψ === nothing
            ψ = x -> SA[0.0, 0.0, 0.0]
        end
        return HexPattern(α, qs, phases, ψ)
    elseif N == 2
        if ψ === nothing
            ψ = x -> SA[0.0, 0.0]
        end
        return SquarePattern(α, qs, phases, ψ)
    else
        throw(DomainError(N, "N must be 2 or 3"))
    end
end

function gen_splines(xs, ys, hs)
    [Interpolations.CubicSplineInterpolation((xs, ys), h',
                                             extrapolation_bc=Interpolations.Flat())
     for h in hs]
end

function InterpolatedPattern(; xs, ys, z, mag)
    dx = (xs[end] - xs[1]) / length(xs)
    dy = (ys[end] - ys[1]) / length(ys)
    Hₚs = gen_splines(xs, ys, H(mag, z, dx, dy))
    InterpolatedPattern(SA[0.0, 0.0], SVector{3}([Hₚs[i] for i in 1:3]))
end

Gₙ(x, y, z) = @. 1 / ((x^2 + y^2 + z^2)^(3 / 2))
restrict(r, Δr, v=0) = ifelse.(abs.(r) .> Δr / 2, v, r)

function H(mag, z, dx, dy)
    ny, nx = size(mag)
    xs = FFTW.ifftshift(range(-nx / 2 * dx, nx / 2 * dx - dx, length=nx))
    ys = FFTW.ifftshift(range(-ny / 2 * dy, ny / 2 * dy - dy, length=ny))
    Δr = min(ny * dy, nx * dx)
    xs = restrict(xs, Δr, 1e10)
    ys = restrict(ys, Δr, 1e10)
    xx, yy = meshgrid(xs, ys)
    gn = Gₙ(xx, yy, z)
    Gx = FFTW.fft(xx .* gn)
    Gy = FFTW.fft(yy .* gn)
    Gz = FFTW.fft(z .* gn)
    Mz = FFTW.fft(mag)
    Hx = real.(FFTW.ifft(Gx .* Mz))
    Hy = real.(FFTW.ifft(Gy .* Mz))
    Hz = real.(FFTW.ifft(Gz .* Mz))
    Hz = Hz .- mean(Hz)
    scale = 1 / mean(abs.(Hz))
    Hx = Hx .* scale
    Hy = Hy .* scale
    Hz = Hz .* scale
    return Hx, Hy, Hz
end

function M(r, pat::AnalyticPattern)
    qr = pat.qs * r - pat.phases - pat.ψ(r)
    amp = 0.0
    if typeof(pat) <: HexPattern
        amp = cos(3.0 * mean(pat.ψ(r))) / 2.0
    end
    for i in 1:N(pat)
        @inbounds amp += cos(qr[i])
    end
    sign(amp)
end
