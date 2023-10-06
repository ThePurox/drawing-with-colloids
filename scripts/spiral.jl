using DrWatson
@quickactivate

using StaticArrays
using AdaptiveBD
include(srcdir("spiral.jl"))

ψs(r) = SA[0.0, 0.0, atan(r[1], r[2])]
pat = genPattern(0.0, 3, SVector{2}([0.0, 0.0]), ψs)

T() = 10.0
control_space(t) = spiral_loop(t)
xra = -50:0.03:50

mag = [M(SA[x, y], pat) for y in xra, x in xra];
cpat = InterpolatedPattern(mag=mag, xs=xra, ys=xra, z=0.5);
s = System(ps=random_particles(1), T=0.0, γ=1.0, t=0.0, interactions=(),
           externalForces=(fmag!,), pats=[cpat]);
rr = RSwM(s=deepcopy(s), relTol=1e-2, absTol=1e-4);
data_in = solve!(rr, 200 * T(), T() / 200);
data_out = solve!(rr, 400 * T(), T() / 200);

plot_pats(data_in[2], color="tab:blue")
plot_pats(data_out[2], color="tab:orange")
