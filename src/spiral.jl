using DrWatson
@quickactivate

using DataStructures
using Statistics

include("plot-tools.jl")
include("tools.jl")
include("forces.jl")
include("loops.jl")
include("gen-patterns.jl")

function spiral_loop(t)
    nᵢₙ = 200
    if t < nᵢₙ * T()
        return spiral_in(t)
    else
        return spiral_out(t)
    end
end
