module AdaptiveBD
using StaticArrays
using DataStructures
using LinearAlgebra
using ProgressMeter
using Infiltrator

export Particle, BareParticle, System, RSwM, step!, solve!

mutable struct Particle
    r::SVector{2, Float64}
    v::SVector{2, Float64}
    F::SVector{2, Float64}
    rRand::SVector{2, Float64}
    χ::Float64
end

struct BareParticle
    r::SVector{2, Float64}
    χ::Float64
end

BareParticle(p::Particle) = BareParticle(p.r, p.χ)

function Particle(r::AbstractVector, χ::Float64=1.0)
    r = SVector{2}(r)
    v = SVector{2}([0.0, 0.0])
    F = SVector{2}([0.0, 0.0])
    rRand = SVector{2}([0.0, 0.0])
    Particle(r, v, F, rRand, χ)
end

mutable struct System{F, G, P}
    ps::Vector{Particle}
    T::Float64
    γ::Float64
    stdthermal::Float64
    t::Float64
    interactions::F
    externalForces::G
    pats::P
    Hₑ::SVector{3, Float64}
end

function System(; ps::Vector{Particle},
                T::Float64,
                γ::Float64,
                t::Float64,
                interactions::Tuple,
                externalForces::Tuple,
                pats)
    System(ps, T, γ, sqrt(2T / γ), t, interactions, externalForces, pats, SA[0.0, 0.0, 0.0])
end

struct RandomIncrement
    R::Vector{SVector{2, Float64}}
    dt::Float64
end
RandomIncrement(r::RandomIncrement) = RandomIncrement(r.R, r.dt)

mutable struct RSwM{F, G, P}
    s::System{F, G, P}
    dtTrial::Float64
    dt::Float64
    dtCutoff::Float64
    relTol::Float64
    absTol::Float64
    qmin::Float64
    qmax::Float64
    α::Float64
    β::Float64

    error::Float64
    q::Float64
    rejections::Int64

    R::Vector{SVector{2, Float64}}
    Rbridge::Vector{SVector{2, Float64}}
    Rdiff::Vector{SVector{2, Float64}}
    Rtmp::Vector{SVector{2, Float64}}

    rsBefore::Vector{SVector{2, Float64}}
    drs₀::Vector{SVector{2, Float64}}
    drs₁::Vector{SVector{2, Float64}}

    Sf::Stack{RandomIncrement}
    Su::Stack{RandomIncrement}

    rejection_counter::Int64
    accept_counter::Int64
end

function RSwM(; s::System, dtTrial::Float64=1e-5, relTol::Float64=1e-4,
              absTol::Float64=1e-6, qmin::Float64=1e-3, qmax::Float64=1.125, α::Float64=2.0,
              β::Float64=0.9, dtCutoff::Float64=1e-13)
    N = length(s.ps)
    R = [zeros(SVector{2}) for _ in 1:N]
    Rbridge = [zeros(SVector{2}) for _ in 1:N]
    genRandVec!(R, 0.0, dtTrial)
    Rdiff = [zeros(SVector{2}) for _ in 1:N]
    Rtmp = [zeros(SVector{2}) for _ in 1:N]
    rsBefore = [zeros(SVector{2}) for _ in 1:N]
    drs₀ = [zeros(SVector{2}) for _ in 1:N]
    drs₁ = [zeros(SVector{2}) for _ in 1:N]

    Sf = Stack{RandomIncrement}()
    Su = Stack{RandomIncrement}()

    rejections = 0
    error = 0.0
    q = 0.0

    RSwM(s, dtTrial, dtTrial, dtCutoff, relTol, absTol, qmin, qmax, α, β, error, q,
         rejections, R, Rbridge, Rdiff, Rtmp, rsBefore, drs₀, drs₁, Sf, Su,
         0, 0)
end
setZero!(R::Vector{SVector{2, Float64}}) =
    for i in 1:length(R)
        R[i] = @SVector zeros(2)
    end

saveParticlesBefore!(r::RSwM) =
    for (i, p) in enumerate(r.s.ps)
        r.rsBefore[i] = p.r
    end
restoreParticlesBefore!(r::RSwM) =
    for (i, p) in enumerate(r.s.ps)
        p.r = r.rsBefore[i]
    end

function genRandVec!(r::RSwM)
    genRandVec!(r.R)
    nothing
end

function genRandVec!(buf::Vector{SVector{2, Float64}}, mean::Float64=0.0, std::Float64=1.0)
    for i in 1:length(buf)
        buf[i] = (@SVector randn(2)) * sqrt(std) .+ mean
    end
end

function genRandVecBridge!(buf::Vector{SVector{2, Float64}},
                           inc::Vector{SVector{2, Float64}}, q::Float64, dt::Float64)
    for i in 1:length(buf)
        buf[i] = (@SVector randn(2)) * sqrt(q * (1 - q) * dt) + q * inc[i]
    end
end

function maybePush!(r::RSwM, S::Stack{RandomIncrement}, inc::RandomIncrement)
    if inc.dt > r.dtCutoff
        push!(S, inc)
    end
end

function applyForces!(r::RSwM, dr::Vector{SVector{2, Float64}})
    for (i, p) in enumerate(r.s.ps)
        # p.rRand = sqrt(2r.s.T/r.s.γ) * r.R[i]
        p.rRand = r.s.stdthermal * r.R[i]
        p.v = p.F / r.s.γ
        dr[i] = p.v * r.dtTrial + p.rRand
        p.r += dr[i]
    end
end

function error!(r::RSwM)
    err = 0.0
    for i in 1:length(r.s.ps)
        err = max(norm(r.drs₀[i] - r.drs₁[i]) /
                  (r.absTol + 0.5 * norm(r.drs₀[i] + r.drs₁[i]) * r.relTol), err)
    end
    r.error = err
    nothing
end

function adaptStepSize!(r::RSwM)
    needSmallerStep = false
    error!(r)
    r.q = 1.0 / (r.α * r.error)^2
    if r.q < 1
        needSmallerStep = true
        r.rejections += 1
        rejectStep!(r)
    else
        needSmallerStep = false
        r.dt = r.dtTrial
        prepareNextStep!(r)
    end
    return needSmallerStep
end

function rejectStep!(r::RSwM)
    r.rejection_counter += 1
    r.q = max(r.qmin, r.q)
    dts = 0.0
    setZero!(r.Rtmp)
    while !isempty(r.Su)
        L = first(r.Su)
        dts += L.dt
        r.Rtmp += L.R
        if dts < (1.0 - r.q) * r.dtTrial
            push!(r.Sf, L)
            pop!(r.Su)
        else
            dtu = L.dt
            dtM = dts - (1.0 - r.q) * r.dtTrial
            qM = dtM / dtu
            genRandVecBridge!(r.Rbridge, L.R, qM, dtu)
            r.Rdiff = L.R - r.Rbridge
            r.R += r.Rbridge - r.Rtmp
            maybePush!(r, r.Sf, RandomIncrement(r.Rdiff, (1.0 - qM) * dtu))
            pop!(r.Su)
            push!(r.Su, RandomIncrement(r.Rbridge, qM * dtu))
            break
        end
    end
    r.dtTrial *= r.q
    nothing
end

function prepareNextStep!(r::RSwM)
    r.accept_counter += 1
    r.dtTrial *= min(r.qmax, r.β * r.q)
    empty!(r.Su)
    dts = 0.0
    setZero!(r.R)
    while !isempty(r.Sf)
        L = first(r.Sf)
        dtf = L.dt
        if dts + dtf <= r.dtTrial
            dts += dtf
            r.R += L.R
            push!(r.Su, L)
            pop!(r.Sf)
        else
            qM = (r.dtTrial - dts) / dtf
            genRandVecBridge!(r.Rbridge, L.R, qM, dtf)
            r.Rdiff = L.R - r.Rbridge
            r.R += r.Rbridge
            pop!(r.Sf)
            maybePush!(r, r.Sf, RandomIncrement(r.Rdiff, (1.0 - qM) * dtf))
            dts += qM * dtf
            break
        end
    end
    if r.dtTrial - dts > r.dtCutoff
        genRandVec!(r.Rtmp, 0.0, r.dtTrial - dts)
        r.R .+= r.Rtmp
        push!(r.Su, RandomIncrement(r.Rtmp, r.dtTrial - dts))
    end
end

function clearForces!(p::Particle)
    p.F = @SArray [0.0, 0.0]
    p.v = @SArray [0.0, 0.0]
    return nothing
end

interact!(s::System) =
    for int in s.interactions
        for i in 1:length(s.ps)
            for j in 1:(i - 1)
                int(s.ps[i], s.ps[j])
            end
        end
    end

function updateForces!(s::System)
    for p in s.ps
        clearForces!(p)
    end
    for ext in s.externalForces
        ext(s)
    end
    interact!(s)
end

function step!(r::RSwM)
    r.rejections = 0
    saveParticlesBefore!(r)
    # updateH!(r.s)

    while true
        if r.rejections > 0
            restoreParticlesBefore!(r)
        end

        updateForces!(r.s)
        for i in 1:length(r.s.ps)
            r.s.ps[i].rRand = r.s.stdthermal * r.R[i]
            r.s.ps[i].v = r.s.ps[i].F / r.s.γ
            r.drs₀[i] = r.s.ps[i].v * r.dtTrial + r.s.ps[i].rRand
            r.s.ps[i].r += r.drs₀[i]
        end
        updateForces!(r.s)
        for i in 1:length(r.s.ps)
            r.s.ps[i].v = r.s.ps[i].F / r.s.γ
            r.drs₁[i] = r.s.ps[i].v * r.dtTrial + r.s.ps[i].rRand
        end
        adaptStepSize!(r) || break
    end
    for i in 1:length(r.s.ps)
        r.s.ps[i].r += (r.drs₁[i] - r.drs₀[i]) / 2
    end
    r.s.t += r.dt
    return nothing
end

function solve!(r::RSwM, t_final::Float64, sample_every::Float64)
    sample_times = (r.s.t):sample_every:(t_final + sample_every)
    p = Progress(length(sample_times))
    hist = Queue{Float64}()
    time = Queue{Float64}()
    parts = Queue{Vector{BareParticle}}()

    next_sample, state = iterate(sample_times, 0)

    try
        while r.s.t < t_final
            try
                if r.s.t > next_sample
                    enqueue!(hist, r.dt)
                    enqueue!(parts, [BareParticle(p) for p in r.s.ps])
                    enqueue!(time, r.s.t)
                    next_sample, state = iterate(sample_times, state)
                    next!(p)
                end
            catch
                @infiltrate
                println("could not sample")
                println("t: $(r.s.t)")
                println("ps: $(r.s.ps)")
            end

            step!(r)
        end
    catch
        println("simulation did not finish!")
        println("an error occured at time $(r.s.t)")
    end
    enqueue!(hist, r.dt)
    enqueue!(parts, [BareParticle(p) for p in r.s.ps])
    enqueue!(time, r.s.t)
    finish!(p)
    hist, parts, time
end

end
