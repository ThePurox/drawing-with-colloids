using DrWatson
@quickactivate

T() = 5.0
const Nₚ = 18
const dϕ = 1 / (2Nₚ)
const dϕi = 0.0275
const height = 20

include(srcdir("alphabet.jl"))

ps = [genPattern(π / 4 + i * π * dϕ, 2, SA[0.0, 0.0]) for i in 0:(Nₚ - 1)]

m_funcs = [(x, y) -> sign(M(x, y, p)) for p in ps]

xs = -50:0.03:50
ys = -150.02:0.03:-50.03
xx, yy = meshgrid(xs, ys)

println("Generating Patterns, this will take some minutes")
ms = gen_mags()
cpats = gen_cpats(ms)

control_space(t) = write_alphabet_compact(t)
ss = [System(ps=[Particle([0.0, -100.0])], T=0.0, γ=1.0, t=0.0, interactions=(),
             externalForces=(fmag!,), pats=[cpat]) for cpat in cpats];
rrs = [RSwM(s=s, dtTrial=1e-4, relTol=1e-2, absTol=1e-4) for s in ss];
datas = Folds.map(r -> solve!(r, length(cc) * T(), T() / 100), rrs);
plot_all_letters(datas)
