using DrWatson
@quickactivate
using Random

include(srcdir("ps-symbols.jl"))

xs = -50:0.03:50
ys = -150.02:0.03:-50.03
ystotal = -150.02:0.03:50

ψs(r, pat::HexPattern) = SA[0.0, 0.0, atan(r[1], r[2])]
ry(ψs) = √3 * (ψs[3] / 3 - (ψs[1] + ψs[2]) / 6)

T() = 5.0

control_space(t) = ps_loop(t)

println("Generating patterns, this may take some minutes")
cpats = [InterpolatedPattern(mag=mag, xs=xs, ys=ystotal, z=0.5) for mag in gen_mags()]

function sim_symbols(Tl=0.0)
    ss = [System(ps=random_particles(1), T=Tl, γ=1.0, t=0.0, interactions=(),
                 externalForces=(fmag!,), pats=[c]) for c in cpats]
    rrs = [RSwM(s=deepcopy(s), absTol=1e-4, relTol=1e-2) for s in ss]
    data = Folds.map(x -> solve!(x, 350.0T(), T() / 200.0), rrs)
    fig, axs = subplots(1, 3)
    for (d, ax) in zip(data, axs)
        plot_pats(d[2]; ax)
    end
    data
end

data = sim_symbols();
