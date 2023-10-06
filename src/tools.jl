using AdaptiveBD
using LinearAlgebra
using StaticArrays

random_particle(r=30.0) = Particle(SA[-r, -r] + 2r * rand(2))
random_particles(N, r=30.0) = [random_particle(r) for _ in 1:N]

function updateH!(s::System, c)
    ϕ, θ = c(s.t)
    s.Hₑ = SA[cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)]
    nothing
end

function meshgrid(xs, ys)
    xx = xs' .* ones(size(ys)[1])
    yy = ones(size(xs)[1])' .* ys
    return xx, yy
end
