using DrWatson
@quickactivate

using Test
include(srcdir("loops.jl"))

T() = 5.0

@testset "loops" begin
    f(t) = loopN(4, SA[(0.0, 0.0), (1, π)], t)
    @test all(f(0.0) == f(T()) .≈ (2π / 8, 0.0))
    @test all(f(1T() / 8) .≈ (2π * 2 / 8, 0.0))
    @test all(f(T() / 4) .≈ (2π * 3 / 8, 0.0))
    @test all(f(3T() / 8) .≈ (2π * 3 / 8, π / 2))
    @test all(f(T() / 2) .≈ (2π * 3 / 8, π))
    @test all(f(5T() / 8) .≈ (2π * 2 / 8, π))
    @test all(f(3T() / 4) .≈ (2π / 8, π))
    @test all(f(7T() / 8) .≈ (2π / 8, π / 2))

    @testset "loop4p" for t in [i * T() / 8 for i in 1:8]
        @test loop4p(0.5, SA[(0.0, 0.0), (1.0, π)], t) == loop4(SA[(0.5, 0.0), (1.5, π)], t)
    end
end
