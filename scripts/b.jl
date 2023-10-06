using DrWatson
@quickactivate

include(srcdir("fromimage.jl"))
include(srcdir("forces.jl"))
include(srcdir("plot-tools.jl"))
include(srcdir("loops.jl"))

T() = 8.0
control_loop(t) = b_loop(t)
xra = -50:0.03:50
xras = -50:50

println("generating pattern")
raw_image = load(projectdir("assets/b.png"))
image = HSVA.(raw_image[end:-1:1, :])
rs, tree = gentree(image)
pc, ps = gen_interpolated_phase(rs, image);
pat = genPattern(0.0, 3, SVector{2}([0.0, 0.0]), ψs)
mag, cpat, s = mcs();

println("simulating B")
Np = 40
rr = RSwM(s=deepcopy(s), relTol=1e-2, absTol=1e-4);
data = solve!(rr, 600T(), T() / 200);
fig, ax = subplots()
ax.imshow(mag, origin="lower", extent=[extrema(xra)..., extrema(xra)...], cmap="binary")
plot_pats(data[2], ax=ax)
