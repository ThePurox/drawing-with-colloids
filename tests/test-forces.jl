using DrWatson
@quickactivate

using Test
include(srcdir("forces.jl"))

control_space(t) = SA[π / 2 - π / 2 * cos(2π * t), π / 2 * sin(2π * t)]
@testset "Magnetic Forces" begin
    sq_pat = genPattern(0.0, 2, SA[0.0, 0.0])
    s = System(ps=[Particle([0.0, 0.0])], T=0.0, γ=1.0, t=0.0, interactions=(),
               externalForces=(fmag!,), pats=[sq_pat])
    @testset "fmagz!" begin
        fmag!(s)
        @test s.ps[1].F≈[0.0, 0.0] atol=1e-13
        AdaptiveBD.clearForces!(s.ps[1])
        s.ps[1].r = SA[0.1, 0.1]
        fmag!(s)
        @test all(s.ps[1].F .< 0)
        AdaptiveBD.clearForces!(s.ps[1])
        @test s.ps[1].F≈[0.0, 0.0] atol=1e-13
        s.ps[1].r = SA[-0.1, -0.1]
        fmag!(s)
        @test all(s.ps[1].F .> 0)
        AdaptiveBD.clearForces!(s.ps[1])
    end
    @testset "fmag!" begin
        s.t = 1e-6
        s.ps[1].r = SA[0.0, 0.0]
        fmag!(s)
        @test s.ps[1].F≈[0.0, 0.0] atol=1e-3
        AdaptiveBD.clearForces!(s.ps[1])
        s.ps[1].r = SA[0.1, 0.1]
        fmag!(s)
        @test all(s.ps[1].F .< 0)
        AdaptiveBD.clearForces!(s.ps[1])
        @test s.ps[1].F≈[0.0, 0.0] atol=1e-13
        s.ps[1].r = SA[-0.1, -0.1]
        fmag!(s)
        @test all(s.ps[1].F .> 0)
        AdaptiveBD.clearForces!(s.ps[1])
    end
end
