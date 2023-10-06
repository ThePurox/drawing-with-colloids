using StaticArrays

line_segment(x0, x1, t) = x0 + (x1 - x0) * t
line_segment12(x0, x1, t) = (0.5 + x0 + (x1 - x0) * t) * π / 6
line_segmentN(N, x0, x1, t) = (0.5 + x0 + (x1 - x0) * t) * 2π / N

loop6(loop, t) = loopN(12, loop, t)

reverse_loop(f, t) = f(T() - (t % T()))

function loop4p(p, loop, t)
    loop = map(x -> (x[1] + p, x[2]), loop)
    loop4(loop, t)
end

function loopN(N, loop, t)
    l = length(loop)
    tl = t % T()
    tl *= l / T()
    s = floor(Int64, tl) + 1 # current path segment
    tl %= 1.0
    if tl < 0.5
        return line_segmentN(N, loop[s][1], loop[s % l + 1][1], 2 * tl), loop[s][2]
    else
        return (loop[s % l + 1][1] + 0.5) * 2π / N,
               line_segment(loop[s][2], loop[s % l + 1][2], 2 * tl - 1)
    end
end

function new_spiral(t)
    n = 0.01π
    s = 0.99π
    tt = 3.07 / 8
    ttn = 3.1π / 8
    tn = tt * π
    ts = (1.0 - tt) * π
    tsn = (1.0 - 3.1 / 8) * π
    loops = SA[(0.0, n), (6.0, tsn), (9.0, tn), (11.0, s), (7.0, ts), (3.0, ttn)] # inwards
    loop6(loops, t)
end

function spiral(t, tt, d)
    n = 0.01π
    s = 0.99π
    tn = tt * π
    ts = (1.0 - tt) * π
    if d == 'i'
        loops = SA[(0.0, n), (6.0, ts), (9.0, tn), (11.0, s), (7.0, ts), (3.0, tn)] # inwards
    elseif d == 'o'
        loops = SA[(1.0, n), (5.0, ts), (9.0, tn), (11.0, s), (7.0, ts), (3.0, tn)] # spiral outwards
    elseif d == '6'
        loops = SA[(1.0, n), (5.0, tn), (4.7, ts), (8.48, tn), (10.7, ts), (11.0, s),
                   (7.0, ts), (2.52, tn)] # all 6 directions
    end
    loop6(loops, t)
end

function loop3(loop, t)
    # 1 2 u
    # 2 3 ro
    # 3 1 lo
    loopN(3, [(l1 - 0.5, l2) for (l1, l2) in loop], t) #correct for rotation of pattern
end

function loop4(loop, t)
    # 1 2 -> lu
    # 2 3 -> ru
    # 3 4 -> ro
    # 4 1 -> lo
    # loopN(4, [(l1 - 0.5, l2) for (l1, l2) in loop], t) #correct for rotation of pattern
    loopN(4, map(x -> (x[1] - 0.5, x[2]), loop), t) #correct for rotation of pattern
end

spiral_in(t) = new_spiral(t)
spiral_out(t) = spiral(t, 3.5 / 8, 'o')

const neut6 = (-0.5, 0.425π)
const neut4 = (0.0, 0.425π)
const neut3 = (0.0, 0.425π)
