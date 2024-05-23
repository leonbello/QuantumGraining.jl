using Symbolics
using SymbolicUtils
using QuantumGraining

"""
    struct Correction

A struct representing a correction term.

# Fields
- `prefac`: The prefactor of the correction term.
- `exponent`: The exponent of the correction term.
- `poly`: A vector of coefficients for the polynomial term.
- `order`: The order of the polynomial term.

# Constructors
- `Correction(prefac, exponent, poly=Num[1,], order=length(poly))`: Constructs a `Correction` object with the given parameters. If `poly` is not provided, it defaults to `[1]`. If `order` is not provided, it defaults to the length of `poly`.

"""
struct Correction
    prefac
    exponent
    poly::Vector{Num}
    order::Int64

    function Correction(prefac, exponent, poly=Num[1,], order=length(poly))
        if !(poly isa Array)
            poly = [poly]
        end
        if !isequal(poly[1], 0)
            prefac_norm = poly[1]       # Make sure the polynomial is always normalized such that 1 + aτ + bτ^2 + ...
            poly = [1, poly[2:end]./prefac_norm...]
        end
        new(prefac*prefac_norm, exponent, poly, order)
    end
end

Base.show(io::IO, c::Correction) = begin
    @variables τ 
    print(io, to_symbol(c, τ))
end

"""
    struct ContractionCoefficient

A struct representing a contraction coefficient.

# Fields
- `corrections::Vector{Correction}`: Vector of `Correction` objects.
- `exponents::Vector{Number}`: Vector of exponents.
- `prefacs::Vector{Number}`: Vector of prefactors.
- `polys::Vector{Vector{Number}}`: Vector of polynomials.

# Constructors
- `ContractionCoefficient(exponents, prefacs, polys)`: Constructs a `ContractionCoefficient` object with given exponents, prefactors, and polynomials.
- `ContractionCoefficient(exponents, prefacs)`: Constructs a `ContractionCoefficient` object with given exponents and prefactors. Polynomials are set to `[1]` by default.

"""
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
            error("exponents and prefacs should have the same length!")
        end
        polys = fill([1], size(exponents))
        corrections = [Correction(exponents[i], prefacs[i], polys[i]) for i in eachindex(exponents)]
        new(corrections, exponents, prefacs, polys) 
    end
end

"""
    simplify_contraction(c::ContractionCoefficient)

Simplifies the given `ContractionCoefficient` by applying simplification operations to its components.

# Arguments
- `c::ContractionCoefficient`: The `ContractionCoefficient` to be simplified.

# Returns
A new `ContractionCoefficient` object with simplified components.

"""
function simplify_contraction(c::ContractionCoefficient) 
    simplified_polys = [simplify.(poly; simplify_fractions=false) for poly in c.polys]
    simplified_prefacs = simplify.(c.prefacs; simplify_fractions=false)
    simplified_exponents = simplify.(c.exponents; simplify_fractions=false)
    
    try
        simplified_polys = [simplify.(poly; simplify_fractions=true) for poly in c.polys]
    catch
        simplified_polys = [simplify.(poly; simplify_fractions=false) for poly in c.polys]
    end

    try
        simplified_prefacs = simplify.(c.prefacs; simplify_fractions=true)
    catch
        simplified_prefacs = simplify.(c.prefacs; simplify_fractions=false)
    end
    
    try
        simplified_exponents = simplify.(c.exponents; simplify_fractions=true)
    catch
        simplified_exponents = simplify.(c.exponents; simplify_fractions=false)
    end

    return ContractionCoefficient(simplified_exponents, simplified_prefacs, simplified_polys)
end

"""
    contraction_coeff(left::Int, right::Int, freqs::Array)

Calculates the coefficient of a whole contraction, given the contraction and input frequencies. Calculates equation (7) in the paper.

# Arguments
- `left`: the left-order of the contraction
- `right`: the right-order of the contraction
- `freqs`: array of frequencies to put in each mode

# Returns
- `c`: a contraction coefficient struct, symbolic expression for the contraction coeffeicient.
"""
# Possibly, this should be a `ContractionCoefficient` constructor.
function contraction_coeff(left::Int, right::Int, freqs::Array)
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


"""
    extend_correction(c::Correction, poly::Vector{<:Number})

Extend a correction by normalizing and modifying its polynomial representation.

# Arguments
- `c::Correction`: The correction to be extended.
- `poly::Vector{<:Number}`: The polynomial representation of the correction.

# Returns
- `Correction`: The extended correction.

"""
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
"""
    ==(c1::ContractionCoefficient, c2::ContractionCoefficient)

Check if two `ContractionCoefficient` objects are equal by comparing their `corrections` arrays.

# Arguments
- `c1::ContractionCoefficient`: The first `ContractionCoefficient` object.
- `c2::ContractionCoefficient`: The second `ContractionCoefficient` object.

# Returns
- `true` if the `corrections` arrays of `c1` and `c2` are equal, `false` otherwise.
"""
function ==(c1::ContractionCoefficient, c2::ContractionCoefficient)
    return issetequal(c1.corrections, c2.corrections)
end

"""
    ==(c1::Correction, c2::Correction)

Check if two `Correction` objects are equal by comparing their `order`, `exponent`, and `poly` properties.

# Arguments
- `c1::Correction`: The first `Correction` object.
- `c2::Correction`: The second `Correction` object.

# Returns
- `true` if the `order`, `exponent`, and `poly` properties of `c1` and `c2` are equal, `false` otherwise.
"""
function ==(c1::Correction, c2::Correction)
    if c1.order != c2.order
        return false
    else
        exp_eq = isequal(simplify.(expand(c1.exponent) - expand(c2.exponent)), 0)
        poly_eq = isequal(simplify.(c1.prefac*c1.poly - c2.prefac*c2.poly), 0)
    end
    return exp_eq && poly_eq
end

import Base: +
"""
    +(c1::Correction, c2::Correction)

Addition operator for `Correction` objects.

This function adds two `Correction` objects `c1` and `c2` and returns a new `ContractionCoefficient` object.

# Arguments
- `c1::Correction`: The first `Correction` object.
- `c2::Correction`: The second `Correction` object.

# Returns
- `ContractionCoefficient`: The result of adding `c1` and `c2`.

"""
function  +(c1::Correction, c2::Correction)
        if isequal(expand(c1.exponent), expand(c2.exponent))
            exponents = [c1.exponent]
            polys = ordered_sum(c1.prefac.*c1.poly, c2.prefac.*c2.poly)
            if !isequal(polys[1], 0)
                try
                    prefacs = [simplify.(polys[1]; simplify_fractions=true)]
                catch
                    prefacs = [simplify.(polys[1]; simplify_fractions=false)]
                end
                try
                    polys = simplify.(polys/polys[1]; simplify_fractions=true)
                catch
                    polys = simplify.(polys/polys[1]; simplify_fractions=false)
                end
            end
        else
            exponents = [c1.exponent, c2.exponent]
            prefacs = [c1.prefac, c2.prefac]
            polys = [c1.poly, c2.poly]
        end
    return ContractionCoefficient(exponents, prefacs, polys)
end

"""
    +(c1::ContractionCoefficient, c2::Correction)

Add a `Correction` object to a `ContractionCoefficient` object.

# Arguments
- `c1::ContractionCoefficient`: The `ContractionCoefficient` object to add to.
- `c2::Correction`: The `Correction` object to add.

# Returns
A new `ContractionCoefficient` object with the `Correction` object added.

"""
function +(c1::ContractionCoefficient, c2::Correction)
    new_c = deepcopy(c1)
    prefacs, polys, exponents = new_c.prefacs, new_c.polys, new_c.exponents
    push!(exponents, c2.exponent)
    push!(prefacs, c2.prefac)
    push!(polys, c2.poly)
    return merge_duplicate_exponents(new_c)
end

"""
    +(c1::ContractionCoefficient, c2::ContractionCoefficient)

Addition operator for `ContractionCoefficient` objects.

This function performs element-wise addition of two `ContractionCoefficient` objects, `c1` and `c2`. It creates a new `ContractionCoefficient` object by adding the corresponding elements of `c1` and `c2`. 
The resulting object is then simplified by merging duplicate exponents and simplifying the contraction.

# Arguments
- `c1::ContractionCoefficient`: The first `ContractionCoefficient` object.
- `c2::ContractionCoefficient`: The second `ContractionCoefficient` object.

# Returns
- `ContractionCoefficient`: The result of the addition operation.

"""
function  +(c1::ContractionCoefficient, c2::ContractionCoefficient)
    new_c = deepcopy(c1)
    for i in eachindex(c2.exponents)
        corr = Correction(c2.prefacs[i], c2.exponents[i], c2.polys[i])
        new_c += corr
    end
    return simplify_contraction(merge_duplicate_exponents(new_c))
end

import Base: *
"""
    *(c1::Correction, c2::Correction)

Multiply two `Correction` objects.

# Arguments
- `c1::Correction`: The first `Correction` object.
- `c2::Correction`: The second `Correction` object.

# Returns
- `Correction`: The result of multiplying `c1` and `c2`.

"""
function *(c1::Correction, c2::Correction)
    exponent = c1.exponent + c2.exponent
    poly = conv(c1.poly, c2.poly)
    order = c1.order + c2.order - 1
    prefac = c1.prefac*c2.prefac
    return Correction(prefac, exponent, poly, order)
end

"""
    *(c::Correction, n::Number)

Multiply a `Correction` object `c` by a number `n`.

# Arguments
- `c::Correction`: The `Correction` object to be multiplied.
- `n::Number`: The number to multiply `c` by.

# Returns
- `Correction`: The result of multiplying `c` by `n`.

"""
function *(c::Correction, n::Number)
    exponent = c.exponent
    poly = c.poly
    order = c.order
    prefac = c.prefac*n
    return Correction(prefac, exponent, poly, order)
end

"""
    *(n::Number, c::Correction)

Multiply a number `n` with a `Correction` object `c`.

# Arguments
- `n::Number`: The number to be multiplied.
- `c::Correction`: The `Correction` object to be multiplied.

# Returns
- `result`: The result of multiplying `n` with `c`.

"""
*(n::Number, c::Correction) = c*n

"""
    *(n::Number, c1::ContractionCoefficient)

Multiply a number `n` with a `ContractionCoefficient` object `c1`.

# Arguments
- `n::Number`: The number to be multiplied.
- `c1::ContractionCoefficient`: The `ContractionCoefficient` object to be multiplied.

# Returns
- `result`: The result of multiplying `n` with `c1`.

"""
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

"""
    conv(u::Vector{<:Number}, v::Vector{<:Number})

Compute the convolution of two vectors `u` and `v`.

# Arguments
- `u::Vector{<:Number}`: The first input vector.
- `v::Vector{<:Number}`: The second input vector.

# Returns
- `result::Vector{<:Number}`: The resulting vector after performing the convolution.

"""
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

"""
    pad(a::Vector, new_length)

Pad a vector `a` with zeros to a new length `new_length`.

# Arguments
- `a::Vector`: The input vector to be padded.
- `new_length`: The desired length of the padded vector.

# Returns
- A new vector with length `new_length`, where the original elements of `a` are preserved and the remaining elements are filled with zeros.

"""
function pad(a::Vector, new_length)
    old_length = length(a)
    if new_length < old_length
        error("New length must be bigger than old length!")
    end
    
    padding = zeros(typeof(a[1]), new_length - old_length)
    return [a..., padding...]
end

"""
    ordered_sum(a, b)

Compute the element-wise sum of two arrays `a` and `b`, padding the shorter array with zeros.

# Arguments
- `a`: The first input array.
- `b`: The second input array.

# Returns
The element-wise sum of `a` and `b`, with the shorter array padded with zeros.
"""
function ordered_sum(a, b)
    new_length = maximum(length.([a, b]))
    padded_a = pad(a, new_length)
    padded_b = pad(b, new_length)

    return padded_a .+ padded_b
end

"""
    merge_duplicate_exponents(c::ContractionCoefficient)

Merge duplicate exponents in a `ContractionCoefficient` object.

# Arguments
- `c::ContractionCoefficient`: The `ContractionCoefficient` object to merge duplicate exponents.

# Returns
A new `ContractionCoefficient` object with merged duplicate exponents.

# Description
This function takes a `ContractionCoefficient` object and merges any duplicate exponents. It iterates over the exponents, prefactors, and polynomials in the `ContractionCoefficient` object and checks for duplicates. 
If duplicates are found, the prefactors and polynomials are merged into a single polynomial. 
The resulting unique exponents, merged prefactors, and merged polynomials are then used to create a new `ContractionCoefficient` object.

Note that the function uses the `ordered_sum` function to merge the polynomials and the `simplify` function to simplify the resulting polynomial.

"""
function merge_duplicate_exponents(c::ContractionCoefficient)
    unique_exponents = []
    unique_prefacs = []
    unique_polys = []
    merged_indices = []

    for i in eachindex(c.exponents)
        if i in merged_indices
            continue            # skip index if it has already been merged
        end

        exponent = c.exponents[i]
        prefac = c.prefacs[i]
        poly = c.polys[i]

        indices_to_merge = [i]
        #indices_to_merge = []
        for j in (i+1):length(c.exponents)
            if isequal(expand(c.exponents[j]), expand(exponent))
                push!(indices_to_merge, j)
                push!(merged_indices, j)
            end
        end

        if length(indices_to_merge) > 1
            # Merge prefacs and polys for indices_to_merge
            merged_poly = [0]
            for j1 in indices_to_merge
                merged_poly = ordered_sum(merged_poly, c.prefacs[j1]*c.polys[j1])
            end

            try
                merged_poly = simplify.(merged_poly; simplify_fractions=true)
            catch
                merged_poly = simplify.(merged_poly; simplify_fractions=false)
            end

            if !isequal(merged_poly[1], 0)
                merged_prefac = merged_poly[1]
                normalized_merged_poly = [convert(Num, 1)]
                for j2 in 2:length(merged_poly)
                    try
                        global poly_coeff = simplify.(merged_poly[j2]/merged_prefac; simplify_fractions=true)
                    catch
                        global poly_coeff = simplify.(merged_poly[j2]/merged_prefac; simplify_fractions=false)
                    end
                    push!(normalized_merged_poly, poly_coeff)
                end

                push!(unique_exponents, exponent)
                push!(unique_prefacs, merged_prefac)
                push!(unique_polys, normalized_merged_poly)
            end

        else
            push!(unique_exponents, exponent)
            push!(unique_prefacs, prefac)
            push!(unique_polys, poly)
        end
    end

    return ContractionCoefficient(unique_exponents, unique_prefacs, unique_polys)
end
