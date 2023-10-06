
const neut = (0.0, 0.0)
function write_A(α, t)
    h2 = round(Int, height / 2)
    if t / T() < 4h2
        rd(α, t)
    elseif t / T() < 6h2
        lu(α, t)
    elseif t / T() < 8h2
        l(α, t)
    elseif t / T() < 10h2
        ld(α, t)
    elseif t / T() < 14h2 + 1
        ru(α, t)
    else
        neut
    end
end

function write_B(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 1h4
        r(α, t)
    elseif t / T() < 3h4
        rd(α, t)
    elseif t / T() < 5h4
        ld(α, t)
    elseif t / T() < 6h4
        l(α, t)
    elseif t / T() < 7h4
        r(α, t)
    elseif t / T() < 9h4
        rd(α, t)
    elseif t / T() < 11h4
        ld(α, t)
    elseif t / T() < 12h4
        l(α, t)
    elseif t / T() < 16h4
        u(α, t)
    else
        neut
    end
end

function write_C(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        ld(α, t)
    elseif t / T() < 6h4
        d(α, t)
    elseif t / T() < 8h4
        rd(α, t)
    elseif t / T() < 10h4
        r(α, t)
    else
        neut
    end
end

function write_D(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        r(α, t)
    elseif t / T() < 4h4
        rd(α, t)
    elseif t / T() < 6h4
        d(α, t)
    elseif t / T() < 8h4
        ld(α, t)
    elseif t / T() < 10h4
        l(α, t)
    elseif t / T() < 14h4
        u(α, t)
    else
        neut
    end
end

function write_E(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        d(α, t)
    elseif t / T() < 5h4
        r(α, t)
    elseif t / T() < 6h4
        l(α, t)
    elseif t / T() < 8h4
        d(α, t)
    elseif t / T() < 10h4
        r(α, t)
    else
        neut
    end
end

function write_F(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        d(α, t)
    elseif t / T() < 5h4
        r(α, t)
    elseif t / T() < 6h4
        l(α, t)
    elseif t / T() < 8h4
        d(α, t)
    else
        neut
    end
end

function write_G(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        ld(α, t)
    elseif t / T() < 6h4
        d(α, t)
    elseif t / T() < 8h4
        rd(α, t)
    elseif t / T() < 10h4
        r(α, t)
    elseif t / T() < 11.5h4
        u(α, t)
    elseif t / T() < 12.5h4
        l(α, t)
    else
        neut
    end
end

function write_H(α, t)
    h2 = round(Int, height / 2)
    if t / T() < height
        d(α, t)
    elseif t / T() < 3h2
        u(α, t)
    elseif t / T() < 4h2
        r(α, t)
    elseif t / T() < 5h2
        d(α, t)
    elseif t / T() < 7h2
        u(α, t)
    else
        neut
    end
end

function write_I(α, t)
    h4 = round(Int, height / 4)
    if t / T() < h4
        l(α, t)
    elseif t / T() < 3h4
        r(α, t)
    elseif t / T() < 4h4
        l(α, t)
    elseif t / T() < 8h4
        d(α, t)
    elseif t / T() < 9h4
        l(α, t)
    elseif t / T() < 11h4
        r(α, t)
    else
        neut
    end
end

function write_J(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        r(α, t)
    elseif t / T() < 5h4
        d(α, t)
    elseif t / T() < 6h4
        ld(α, t)
    elseif t / T() < 7h4
        l(α, t)
    elseif t / T() < 8h4
        lu(α, t)
    else
        neut
    end
end

function write_K(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 4h4
        d(α, t)
    elseif t / T() < 6h4
        u(α, t)
    elseif t / T() < 10h4
        rd(α, t)
    elseif t / T() < 14h4
        lu(α, t)
    elseif t / T() < 18h4
        ru(α, t)
    else
        neut
    end
end

function write_L(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 4h4
        d(α, t)
    elseif t / T() < 6h4
        r(α, t)
    else
        neut
    end
end

function write_M(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 4h4
        d(α, t)
    elseif t / T() < 8h4
        u(α, t)
    elseif t / T() < 12h4
        rd(α, t)
    elseif t / T() < 16h4
        ru(α, t)
    elseif t / T() < 20h4
        d(α, t)
    else
        neut
    end
end

function write_N(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 4h4
        d(α, t)
    elseif t / T() < 8h4
        u(α, t)
    elseif t / T() < 16h4
        rd(α, t)
    elseif t / T() < 20h4
        u(α, t)
    else
        neut
    end
end

function write_O(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        ld(α, t)
    elseif t / T() < 6h4
        d(α, t)
    elseif t / T() < 8h4
        rd(α, t)
    elseif t / T() < 10h4
        r(α, t)
    elseif t / T() < 12h4
        ru(α, t)
    elseif t / T() < 14h4
        u(α, t)
    elseif t / T() < 16h4
        lu(α, t)
    else
        neut
    end
end

function write_P(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 1h4
        r(α, t)
    elseif t / T() < 3h4
        rd(α, t)
    elseif t / T() < 5h4
        ld(α, t)
    elseif t / T() < 6h4
        l(α, t)
    elseif t / T() < 8h4
        u(α, t)
    elseif t / T() < 12h4
        d(α, t)
    else
        neut
    end
end

function write_Q(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 2h4
        l(α, t)
    elseif t / T() < 4h4
        ld(α, t)
    elseif t / T() < 6h4
        d(α, t)
    elseif t / T() < 8h4
        rd(α, t)
    elseif t / T() < 10h4
        r(α, t)
    elseif t / T() < 11h4
        ru(α, t)
    elseif t / T() < 12h4
        rd(α, t)
    elseif t / T() < 14h4
        lu(α, t)
    elseif t / T() < 15h4
        rd(α, t)
    elseif t / T() < 16h4
        ru(α, t)
    elseif t / T() < 18h4
        u(α, t)
    elseif t / T() < 20h4
        lu(α, t)
    else
        neut
    end
end

function write_R(α, t)
    h4 = round(Int, height / 4)
    if t / T() < 1h4
        r(α, t)
    elseif t / T() < 3h4
        rd(α, t)
    elseif t / T() < 5h4
        ld(α, t)
    elseif t / T() < 6h4
        l(α, t)
    elseif t / T() < 8h4
        u(α, t)
    elseif t / T() < 12h4
        d(α, t)
    elseif t / T() < 14h4
        u(α, t)
    elseif t / T() < 18h4
        rd(α, t)
    else
        neut
    end
end
