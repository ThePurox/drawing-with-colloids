using DrWatson
@quickactivate

using Test
include(srcdir("fromimage.jl"))
test_asset_dir(args...) = projectdir("tests", "test-assets", args...)

@testset "generating phase" begin
    img_data = load(test_asset_dir("single-pixel2x2.png"))
    image = img_data[end:-1:1, :]
    image = HSVA.(image)
    rs, tree = gentree(image)
    @test rs == [[1, 1]]
    myfunc(r) = total_avg(r, rs, image)
    αs = [myfunc(SA[x, y]) for x in 1:2, y in 1:2]
    @test all(isapprox.(αs, -π, rtol=0.05))

    img_data = load(test_asset_dir("double-pixel3x2.png"))
    image = img_data[end:-1:1, :]
    image = HSVA.(image)
    rs, tree = gentree(image)
    @test rs == [[1, 1], [2, 3]]
    myfunc_r2(r) = total_avg_r(r, rs, image)
    myfunc2(r) = total_avg(r, rs, image)
    rrs = [myfunc_r2(SA[x, y]) for x in 1:2, y in 1:3]
    αs = [myfunc2(SA[x, y]) for x in 1:2, y in 1:3]
    @test rrs[1, 1] / norm(rrs[1, 1])≈SA[0, -1] atol=0.08
    @test rrs[1, 2]≈SA[0, -1 / 2] atol=0.07
    @test rrs[1, 3]≈SA[0, 3 / 4] atol=0.05
    @test rrs[2, 1]≈SA[0, -3 / 4] atol=0.07
    @test rrs[2, 2]≈SA[0, 1 / 2] atol=0.05
    @test rrs[2, 3] / norm(rrs[2, 3])≈SA[0, 1] atol=0.05
    @test_throws MethodError total_avg(SA[1.0, 2.0], rs, image)
end
