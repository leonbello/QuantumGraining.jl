using Symbolics
#using SymPy
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
    #@variables τ
    SymPy.@syms τ
    sym = c.prefac*exp(-0.5*τ^2*c.exponent)
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]*(τ^(n-1)) for n in 1:c.order])
    return sym
end

function extend_correction(c::Correction, poly::Vector{<:Number})
    non_zero_poly = findall(x -> x != 0, poly)
    if isempty(non_zero_poly)
        poly = [1]
        norm = 1
    else
        norm = poly[non_zero_poly[1]] # find the first non-zero element of poly and normalize the correction by it.
    end
    return Correction(norm*c.prefac,
                    c.exponent,
                    1/norm*poly,
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
    exp_in = 0
    for expon in exponents
        if isequal(c2.exponent, expon)
            exp_in = 1
        end
    end
    if exp_in == 1
        ind = findfirst(item -> isequal(item,c2.exponent), c1.exponents)
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


function  +(c1::ContractionCoefficient, c2::ContractionCoefficient)
    contractCoeff = c1
    for i in 1:length(c2.exponents)
        corr = Correction(c2.prefacs[i], c2.exponents[i], c2.polys[i])
        contractCoeff += corr
    end
    return contractCoeff
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
function *(n::Number, c1::ContractionCoefficient)
    prefacs = c1.prefacs*n
    return ContractionCoefficient(c1.exponents,prefacs,c1.polys)
end

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