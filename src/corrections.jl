using Symbolics
using QuantumGraining

struct Correction
    prefac
    exponent
    poly::Vector{Num}
    order::Int64

    function Correction(prefac, exponent, poly=Num[1,], order=length(poly))
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
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]*(τ^(n-1)) for n in 1:c.order])
    return sym
end

function extend_correction(c::Correction, poly::Vector{<:Number})
    return Correction(poly[1]*c.prefac,
                    c.exponent,
                    1/poly[1]*poly,
                    length(poly)
                    )
end

import Base: +
function  +(c1::Correction, c2::Correction)
        if c1.exponent ≈ c2.exponent
            exponents = [c1.exponent]
            prefacs = [c1.prefac + c2.prefac]
            polys = [(c1.prefacs*c1.poly + c2*prefacsc2.poly)/prefacs[1]]
        else
            exponents = [c1.exponent, c2.exponent]
            prefacs = [c1.prefac, c2.prefac]
            polys = [c1.poly, c2.poly]
        end
    return ContractionCoefficient(exponents, prefacs, polys)
end

function  +(c1::ContractionCoefficient, c2::Correction)
    prefacs,polys,exponents = c1.prefacs, c1.polys, c1.exponents
    if c2.exponent ∈ exponents
        ind = findfirst(item -> item ≈ c2.exponent, c1.exponents)
        prefac_og = prefacs[ind]
        replace!(prefacs, prefacs[ind] => prefacs[ind] + c2.prefac)
        replace!(polys, polys[ind] => (prefac_og*polys[ind] + c2.prefac*c2.poly)/prefacs[ind])
    else
        push!(exponents, c2.exponent)
        push!(prefacs, c2.prefac)
        push!(polys, c2.poly)
    end
    return ContractionCoefficient(exponents, prefacs, polys)
end


# function +(c1::ContractionCoefficient, c2::ContractionCoefficient)
#     prefacs,polys,exponents = c1.prefacs, c1.polys, c1.exponents
#     common_exponents = intersect(c1.exponents, c2.exponents)
#     for expon in common_exponents
#         ind1 = findfirst(item -> item ≈ expon, c1.exponents)
#         ind2 = findfirst(item -> item ≈ expon, c2.exponents)
#         prefac_og = prefacs[ind1]
#         replace!(prefacs, prefacs[ind1] => prefacs[ind1] + c2.prefacs[ind2])
#         replace!(polys, polys[ind1] => (prefac_og*polys[ind1] + c2.prefacs[ind2]*c2.polys[ind2])/prefacs[ind1])
#         deleteat!(c2.exponents, ind2)
#         deleteat!(c2.prefacs, ind2)
#         deleteat!(c2.polys, ind2)
#     end
#     prefacs,polys,exponents = c1.prefacs, c1.polys, c1.exponents
#     append!(exponents, c2.exponents)
#     append!(prefacs, c2.prefacs)
#     append!(polys, c2.polys)
#     return ContractionCoefficient(exponents, prefacs, polys)
# end


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

function conv(u::Vector{<:Number}, v::Vector{<:Number})
    m, n = length(u), length(v)
    result = zeros(eltype(u), m + n - 1)

    for i in 1:m
        for j in 1:n
            result[i+j-1] += u[i] * v[j]
        end
    end

    return result
end