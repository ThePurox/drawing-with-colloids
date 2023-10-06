using AdaptiveBD
include("pattern.jl")

function fmag!(s::System)
    updateH!(s)
    for pat in s.pats
        for p in s.ps
            fmag!(s, pat, p)
        end
    end
    return nothing
end

function updateH!(s::System)
    ϕ, θ = control_space(s.t)
    s.Hₑ = SA[cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)]
    nothing
end

function fmagz!(s::System, pat::AnalyticPattern, p::Particle)
    qr = pat.qs * p.r - pat.phases - pat.ψ(p.r)
    for i in 1:N(pat)
        amp = -q(pat) * sin(qr[i])
        p.F += 2amp * p.χ * s.Hₑ[3] * pat.qs[i, :]
    end
    nothing
end

function fmag!(s::System, pat::AnalyticPattern, p::Particle)
    if s.Hₑ[3] == 1 || s.Hₑ[3] == -1
        fmagz!(s, pat, p)
    else
        qr = pat.qs * p.r - pat.phases - pat.ψ(p.r)
        for i in 1:N(pat)
            @inbounds amp = pat.qs[i, 1] * cos(qr[i]) * s.Hₑ[1] +
                            pat.qs[i, 2] * cos(qr[i]) * s.Hₑ[2] -
                            q(pat) * sin(qr[i]) * s.Hₑ[3]
            @inbounds p.F += 2amp * p.χ * pat.qs[i, :]
        end
    end
    nothing
end

function fmag!(s::System, pat::InterpolatedPattern, p::Particle)
    for (i, H) in enumerate(pat.Hₚs)
        ∇H = Interpolations.gradient(H, p.r...)
        p.F += p.χ .* ∇H .* s.Hₑ[i]
    end
end
