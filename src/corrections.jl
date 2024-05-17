using Symbolics
using SymbolicUtils
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
            error("exponents, prefacs and polys should have the same length!")
        end
        corrections = [Correction(exponents[i], prefacs[i], polys[i]) for i in eachindex(exponents)]
        new(corrections, exponents, prefacs, polys)
    end

    function ContractionCoefficient(exponents, prefacs)
        if !(length(exponents) == length(prefacs))
            error("exponents, prefacs and polys should have the same length!")
        end
        polys = fill([1], size(exponents))
        corrections = [Correction(exponents[i], prefacs[i], polys[i]) for i in eachindex(exponents)]
        new(corrections, exponents, prefacs, polys)
    end
end

function simplify_contraction(c::ContractionCoefficient)
    simplified_prefacs = simplify.(c.prefacs; simplify_fractions=false)
    simplified_exponents = simplify.(c.exponents; simplify_fractions=false)
    simplified_polys = [simplify.(poly; simplify_fractions=false) for poly in c.polys]

    return ContractionCoefficient(simplified_exponents, simplified_prefacs, simplified_polys)
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
    c = ContractionCoefficient(exp_list, pre_list, poly_list)
    return merge_duplicate_exponents(c)
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
    sym = expand(simplify(sym))
    return sym
end

function to_symbol(c::Correction, τ) #where {T <: Number}
    sym = c.prefac*exp(-0.5*τ^2*c.exponent)
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]*(τ^(n-1)) for n in 1:c.order])
    sym = expand(simplify(sym))
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
            max_length = maximum([length(c1.prefacs*c1.poly), length(c2*prefacsc2.poly)])
            polys = [(vcat(c1.prefacs*c1.poly, fill(0, max_length - length(c1.prefacs*c1.poly))) + vcat(c2*prefacsc2.poly, fill(0, max_length - length(c2*prefacsc2.poly))))/prefacs[1]]
        else
            exponents = [c1.exponent, c2.exponent]
            prefacs = [c1.prefac, c2.prefac]
            polys = [c1.poly, c2.poly]
        end
    return ContractionCoefficient(exponents, prefacs, polys)
end

function  +(c1::ContractionCoefficient, c2::Correction)
    new_c = deepcopy(merge_duplicate_exponents(c1))
    C1Exponents = []
    for i in eachindex(c1.exponents)
        push!(C1Exponents, expand(c1.exponents[i]))
    end

    prefacs, polys, exponents = new_c.prefacs, new_c.polys, new_c.exponents
    found = false
    for exponent in exponents
        if isequal(expand(c2.exponent), expand(exponent))
        if isequal(expand(c2.exponent), expand(exponent))
            found = true
            ind = findfirst(item -> isequal(item, expand(c2.exponent)), expand.(c1.exponents))
            old_prefac = prefacs[ind]
            replace!(prefacs, prefacs[ind] => simplify(prefacs[ind] + c2.prefac, simplify_fractions=false))
            replace!(polys, polys[ind] => (ordered_sum(old_prefac*polys[ind], c2.prefac*c2.poly))/prefacs[ind])
        end
    end
    if !found
        push!(exponents, c2.exponent)
        push!(prefacs, c2.prefac)
        push!(polys, c2.poly)
    end
    return merge_duplicate_exponents(ContractionCoefficient(exponents, prefacs, polys))
end


function  +(c1::ContractionCoefficient, c2::ContractionCoefficient)
    new_c = deepcopy(c1)
    for i in eachindex(c2.exponents)
        corr = Correction(c2.prefacs[i], c2.exponents[i], c2.polys[i])
        new_c += corr
    end
    return simplify_contraction(merge_duplicate_exponents(new_c))
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

function pad(a::Vector, new_length)
    old_length = length(a)
    if new_length < old_length
        error("New length must be bigger than old length!")
    end
    
    padding = zeros(typeof(a[1]), new_length - old_length)
    return [a..., padding...]
end

function ordered_sum(a, b)
    new_length = maximum(length.([a, b]))
    padded_a = pad(a, new_length)
    padded_b = pad(b, new_length)

    return padded_a .+ padded_b
end

function merge_duplicate_exponents(c::ContractionCoefficient)
    unique_exponents = []
    unique_prefacs = []
    unique_polys = []
    merged_indices = []

    for i in eachindex(c.exponents)
        if i in merged_indices
            continue
        end

        exponent = c.exponents[i]
        prefac = c.prefacs[i]
        poly = c.polys[i]

        indices_to_merge = [i]
        #indices_to_merge = []
        for j in (i+1):length(c.exponents)
            if isequal(expand(c.exponents[j]^2) - expand(exponent^2), 0)
                push!(indices_to_merge, j)
                push!(merged_indices, j)
            end
        end

        if length(indices_to_merge) > 1
            # Merge prefacs and polys for indices_to_merge
            merged_prefac = sum(simplify.(c.prefacs[indices_to_merge]))

            merged_poly = [0]
            for i in indices_to_merge
                merged_poly = ordered_sum(c.polys[i], merged_poly)
            end
            merged_poly = simplify.(merged_poly, simplify_fractions=false)
            
            #merged_poly = simplify(sum(c.polys[indices_to_merge]))


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
