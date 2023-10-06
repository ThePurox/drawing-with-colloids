using NearestNeighbors
using ImageCore
using ColorTypes
using StaticArrays
using ProgressMeter
using Statistics

include("pattern.jl")

function gentree(image)
    rs = Vector{SVector{2, Int64}}()
    for i in eachindex(image)
        if image[i].alpha == 1.0
            push!(rs, SVector(Tuple(CartesianIndices(image)[i])))
        end
    end
    rs, KDTree(rs)
end

phase(image, rs, id) = 0.0 / 0.0
extract_phaseHSV(px) = float(hue(px)) / 360.0

function phase(r, tree, rs, image)
    idx, ds = knn(tree, r, 10)
    p(id) = extract_phaseHSV(image[rs[id]...])
    mean_rot(p.(idx))
end

function plot_phase(; kwargs...)
    p(r) = phase(r, tree, rs, image)
    imshow(p.([(SA[x, y] .+ 50) .* 10 for x in xras, y in xras]), origin="lower",
           extent=[extrema(xras)..., extrema(xras)...]; kwargs...)
end

m(r) = phase(r, tree, rs, image)
ψs(r) = SA[0.0, 0.0, (atan(ps(r...), pc(r...)))]
gen_m() = @showprogress [M(SA[x, y], pat) for x in xra, y in xra]

b_loop(t) = spiral(t, 3.0 / 8, '6') .+ (π / 2, 0.0)
control_space(t) = b_loop(t)

function mcs(Np=40)
    println("generating magnetisation")
    mag = gen_m()
    println("interpolating magnetisation")
    cpat = InterpolatedPattern(mag=mag, xs=xra, ys=xra, z=0.50)
    ps = random_particles(Np, 40.0)
    s = System(ps=ps, T=0.0, γ=1.0, t=0.0, interactions=(), externalForces=(fmag!,),
               pats=[cpat])
    mag, cpat, s
end

function sim(Np=40)
    m, c, s = mcs(Np)
    rr = RSwM(s=s, dtTrial=1e-4, relTol=1e-2, absTol=1e-4)
    data = solve!(rr, 400 * T(), T() / 20)
    m, c, s, data
end

mean_rot(αs) = atan(mean(sin.(2π .* αs)), mean(cos.(2π .* αs))) / 2π

function gen_interpolated_phase(rs, image)
    xs = ys = -50:0.5:50
    p = [total_avg(Int.((SA[x, y] .+ 50) .* 10), rs, image) for x in xs, y in ys]
    ps = sin.(p)
    pc = cos.(p)
    [Interpolations.CubicSplineInterpolation((xs, ys), p) for p in [pc, ps]]
end

total_avg(r::SVector{2, Int64}, rs, image) = atan(total_avg_r(r, rs, image)...)

function total_avg_r(r::SVector{2, Int64}, rs, image)
    r₀ = SA[0.0, 0.0]
    for rᵢ in rs
        α = hue(image[rᵢ...])
        r₀ += SA[sind(α), cosd(α)] ./ (norm(r - rᵢ)^2 + 1e-12)
    end
    r₀
end

total_avg(r) = total_avg(Int.((r .+ 50.0) .* 10.0), rs, image)
total_avg_r(r) = total_avg_r(Int.((r .+ 50.0) .* 10.0), rs, image)
total_avg_r(x, y) = total_avg_r(Int.((SA[x, y] .+ 50.0) .* 10.0), rs, image)

r2img(r) = round.(Int, (r .+ 50) .* 10) .+ SA[1, 1]
u(x, y) = cos(hue(image[r2img(SA[x, y])...]))
v(x, y) = sin(hue(image[r2img(SA[x, y])...]))
unzip(vs) = Tuple(eltype(first(vs))[xyz[j] for xyz in vs] for j in eachindex(first(vs)))

function mean_rot(αs, ds)
    res = SA[0.0, 0.0]
    for (α, d) in zip(αs, ds)
        res += SA[cos(2π * α), sin(2π * α)] ./ (d + 1e-9)
    end
    atan(res[2], res[1]) / 2π
end
