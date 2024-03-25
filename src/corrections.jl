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
Base.show(io::IO, c::Correction) = begin
    @variables τ 
    print(io, to_symbol(c, τ))
end

struct ContractionCoefficient 
    corrections::Vector{Correction}
    exponents::Vector{Number}
    prefacs::Vector{Number}
    polys::Vector{Vector{Number}}

    function ContractionCoefficient(exponents, prefacs, polys)
        if !(length(exponents) == length(prefacs) == length(polys))
            error("exponents, prefacs and polys should have the length!")
        end
        corrections = [Correction(exponents[i], prefacs[i], polys[i]) for i in eachindex(exponents)]
        new(corrections, exponents, prefacs, polys)
    end

    function ContractionCoefficient(exponents, prefacs)
        if !(length(exponents) == length(prefacs))
            error("exponents, prefacs and polys should have the length!")
        end
        polys = fill([1], size(exponents))
        corrections = [Correction(exponents[i], prefacs[i], polys[i]) for i in eachindex(exponents)]
        new(corrections, exponents, prefacs, polys) 
    end
end

# Possibly, this should be a `ContractionCoefficient` constructor.
function contraction_coeff(left::Int, right::Int, freqs::Array)
    """
    contraction_coeff(left::Int, right::Int, freqs::Array)

    Calculates the coefficient of a whole contraction, given the contraction and input frequencies. Calculates equation (7) in the paper.

    Arguments:
    - left: the left-order of the contraction
    - right: the right-order of the contraction
    - freqs: array of frequencies to put in each mode

    Returns:
    - c: a contraction coefficient struct, symbolic expression for the contraction coeffeicient.
    """
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)
    exp_list = []
    pre_list = []
    poly_list = []
    d_list = []
    
    for diagram in diagrams
        reverse!(diagram)                               # reversing since Wentao's order is right-to-left, rather than left-to-right
        ω = split_freqs_into_bubbles(freqs, diagram)
        corr = diagram_correction(ω)
        push!(d_list, diagram)
        push!(exp_list, corr.exponent)
        push!(pre_list, corr.prefac)
        push!(poly_list, corr.poly)
    end
    return ContractionCoefficient(exp_list, pre_list, poly_list)
end
contraction_coeff(order::Tuple{Int, Int}, ω::Array) = contraction_coeff(order[1], order[2], ω)

Base.show(io::IO, coeff::ContractionCoefficient) = begin
    @variables τ
    print(io, to_symbol(coeff, τ))
end

function to_symbol(coeff::ContractionCoefficient, τ) #where {T <: Number}
    sym = 0
    for i in 1:length(coeff.prefacs)
        sym += to_symbol(Correction(coeff.prefacs[i], coeff.exponents[i], coeff.polys[i]), τ)
    end
    return sym
end

function to_symbol(c::Correction, τ) #where {T <: Number}
    sym = c.prefac*exp(-0.5*τ^2*c.exponent)
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]*(τ^(n-1)) for n in 1:c.order])
    return sym
end

function extend_correction(c::Correction, poly::Vector{<:Number})
    #non_zero_poly = findall(x -> x != 0, poly)
    non_zero_poly = findall(x -> !(isequal(x, 0)), poly)
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

import Base: ==
function ==(c1::ContractionCoefficient, c2::ContractionCoefficient)
    return issetequal(c1.corrections, c2.corrections)
end

function ==(c1::Correction, c2::Correction)
    if c1.order != c2.order
        return false
    else
        exp_eq = isequal(simplify(expand(c1.exponent) - expand(c2.exponent)), 0)
        poly_eq = isequal(simplify.(c1.prefac*c1.poly - c2.prefac*c2.poly), 0)
    end
    return exp_eq && poly_eq
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

import Base: -
function  -(c1::Correction, c2::Correction)
    return c1 + -1*c2
end
function  -(c1::ContractionCoefficient, c2::Correction)
    return c1 + -1*c2
end
function  -(c1::ContractionCoefficient, c2::ContractionCoefficient)
    return c1 + -1*c2
end
function  -(c::Correction)
    return -1*c
end
function  -(c::ContractionCoefficient)
    return -1*c
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

function merge_duplicate_exponents(c::ContractionCoefficient)
    unique_exponents = []
    unique_prefacs = []
    unique_polys = []
    merged_indices = []

    for i in 1:length(c.exponents)
        if i in merged_indices
            continue
        end

        exponent = c.exponents[i]
        prefac = c.prefacs[i]
        poly = c.polys[i]

        indices_to_merge = [i]
        for j in (i+1):length(c.exponents)
            if isequal(expand(c.exponents[j]^2) - expand(exponent^2), 0)
                push!(indices_to_merge, j)
                push!(merged_indices, j)
            end
        end

        if length(indices_to_merge) > 1
            # Merge prefacs and polys for indices_to_merge
            merged_prefac = simplify(sum(c.prefacs[idx] for idx in indices_to_merge))
            merged_poly = simplify(sum(c.polys[idx] for idx in indices_to_merge))

            push!(unique_exponents, exponent)
            push!(unique_prefacs, merged_prefac)
            push!(unique_polys, merged_poly)
            
        else
            push!(unique_exponents, exponent)
            push!(unique_prefacs, prefac)
            push!(unique_polys, poly)
        end
    end

    return ContractionCoefficient(unique_exponents, unique_prefacs, unique_polys)
end
