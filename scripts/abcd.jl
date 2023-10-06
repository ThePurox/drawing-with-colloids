using DrWatson
@quickactivate

using StaticArrays

T() = 5.0
const Nₚ = 4
const dϕ = 1 / (2Nₚ)
const dϕi = dϕ
const north = 0.0
const south = Float64(π)
const height = 5

include(srcdir("alphabet.jl"))

function write_abcd(t)
    n_all = floor(Int64, t / (Nₚ * T()))
    t0 = T() * n_all
    ta = t % T() + t0
    n_inner = floor(Int, ((t / T()) % Nₚ))
    if n_inner == 0
        write_A(2n_inner * dϕ, ta)
    elseif n_inner == 1
        write_B(2n_inner * dϕ, ta)
    elseif n_inner == 2
        write_C(2n_inner * dϕ, ta)
    elseif n_inner == 3
        write_D(2n_inner * dϕ, ta)
    else
        return neut
    end
end

function plot_all_letters(datas)
    fig, axs = subplots(2, 2, sharex=true)
    for (i, d) in enumerate(datas)
        axs[i].set_aspect("equal")
        plot_pats(d[2], ax=axs[i])
    end
end

function gen_patterns()
    println("Generating Patterns, this will take some minutes")
    global ms = gen_mags()
    global cpats = gen_cpats(ms)
end

ps = [genPattern(π / 4 + i * π * dϕ, 2, SA[0.0, 0.0]) for i in 0:(Nₚ - 1)]
m_funcs = [(x, y) -> sign(M(x, y, p)) for p in ps]

xs = -50:0.03:50
ys = -50.02:0.03:50.03
xx, yy = meshgrid(xs, ys)

control_space(t) = write_abcd(t)

function simulate()
    ss = [System(ps=[Particle([0.0, 0.0])], T=0.0, γ=1.0, t=0.0, interactions=(),
                 externalForces=(fmag!,), pats=[cpat]) for cpat in cpats]
    rrs = [RSwM(s=s, dtTrial=1e-4, relTol=1e-2, absTol=1e-4) for s in ss]
    datas = Folds.map(r -> solve!(r, 110T(), T() / 100), rrs)
    plot_all_letters(datas)
    datas
end

gen_patterns();
datas = simulate();
