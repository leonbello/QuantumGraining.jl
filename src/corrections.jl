using Symbolics

struct Correction
    prefac
    exponent
    poly::Vector{Float64}
    order::Int64

    function Correction(prefac, exponent, poly=[1,], order=length(poly))
        if !(poly isa Array)
            poly = [poly]
        end
        new(prefac, exponent, poly, order)
    end
end
Base.show(io::IO, c::Correction) = print(io, to_symbol(c))

function to_symbol(c::Correction)
    @variables τ
    sym = c.prefac*exp(-0.5*τ^2*c.exponent)
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]^(τ^(n-1)) for n in 1:c.order])
    return sym
end

import Base: *
function *(c1::Correction, c2::Correction)
    exponent = c1.exponent + c2.exponent
    poly = conv(c1.poly, c2.poly)
    order = c1.order + c2.order - 1
    prefac = c1.prefac*c2.prefac
    return Correction(prefac, exponent, poly, order)
end

function *(c::Correction, n::Number)
    exponent = c.exponent
    poly = c.poly
    order = c.order
    prefac = c.prefac*n
    return Correction(prefac, exponent, poly, order)
end
*(n::Number, c::Correction) = c*n