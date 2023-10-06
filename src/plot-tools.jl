using PyPlot
using Statistics
using LinearAlgebra

include("forces.jl")

function plot_pats(pats; maxparts=typemax(Int), shift=(0.0, 0.0), ax=nothing, lw=3,
                   onlyparts=nothing, skip=0, plotlast=true, kwargs...)
    locpats = hcat([ps for ps in pats]...)
    for i in (onlyparts === nothing ? (1:min(size(locpats)[1], maxparts)) : onlyparts)
        if ax === nothing
            PyPlot.plot([p.r[1] + shift[1] for p in locpats[i, (1 + skip):end]],
                        [p.r[2] + shift[2] for p in locpats[i, (1 + skip):end]], lw=lw;
                        kwargs...)
        else
            ax.plot([p.r[1] + shift[1] for p in locpats[i, (1 + skip):end]],
                    [p.r[2] + shift[2] for p in locpats[i, (1 + skip):end]], lw=lw;
                    kwargs...)
        end
    end
    if plotlast
        for i in (onlyparts === nothing ? (1:min(size(locpats)[1], maxparts)) : onlyparts)
            if ax === nothing
                PyPlot.plot([p.r[1] + shift[1] for p in locpats[i, end:end]],
                            [p.r[2] + shift[2] for p in locpats[i, end:end]], marker="x",
                            markersize=10; kwargs...)
            else
                ax.plot([p.r[1] + shift[1] for p in locpats[i, end:end]],
                        [p.r[2] + shift[2] for p in locpats[i, end:end]], marker="x",
                        markersize=10;
                        kwargs...)
            end
        end
    end
end
