using DrWatson
@quickactivate

using Test
include(srcdir("pattern.jl"))

@testset "pattern" begin
    @testset "square" begin
        @test all(isapprox.(q(1, 2, 0), 2π * SA[-1.0, 0.0], atol=1e-12))
        @test all(isapprox.(q(2, 2, 0), 2π * SA[0.0, -1.0], atol=1e-12))
        @test all(isapprox.(q(1, 2, π / 2), q(2, 2, 0), atol=1e-12))
    end

    @testset "hexagonal" begin
        @test all(isapprox.(q(1, 3, 0), 2π * SA[-1.0, -1 / √3], atol=1e-12))
        @test all(isapprox.(q(2, 3, 0), 2π * SA[1.0, -1 / √3], atol=1e-12))
        @test all(isapprox.(q(3, 3, 0), 4π / √3 * SA[0.0, 1.0], atol=1e-12))
        @test norm(sum(q(i, 3, 0) for i in 1:3))≈0 atol=1e-12
    end
end
