"""
    BVector{T} <: AbstractVector{T}
A struct representing the up- or down-modes. The poles are stored in the poles field.

Arguments:
    - freqs: A vector representing the frequencies of the vector.
    - poles: A vector of the indices of poles of the vector, if any. Does not include the special mode, if there is one.
    - special: A boolean representing whether the vector is the special mode.
    - type: A symbol representing the type of the vector. Must be :up or :down.
"""
mutable struct BVector{T} <: AbstractVector{T}
    freqs::Vector{T}
    poles::Vector{Int}
    special::Bool
    type::Symbol

    function BVector(vec::Vector{T}, type::Symbol, special::Bool) where {T}
        if type != :down && type != :up
            error("type must be :up or :down")
        end
        start = (special == true) ? 2 : 1
        poles = find_poles(vec[start:end])
        return new{T}(vec, poles, special, type)
    end
end

""""
    UVec(u::Vector; special=false)
Creates a BVector with type :up. Assumes the frequency ordering is clockwise
"""
function UVec(u::Vector; special=false)
    return BVector(u, :up, special)
end

""""
    DVec(u::Vector; special=false)
Creates a BVector with type :down. Assumes the frequency ordering is clockwise
"""
function DVec(u::Vector; special=false)
    return BVector(u, :down, special)
end

"""
    vec_factorial(u; include_poles=true)
Calculates vector factorial for a BVector u. If the vector includes the edge mode, then the first term is omitted.
"""
function vec_factorial(u::BVector; include_poles=false)
    if u.type == :down
        reverse!(u.freqs)
    end    

    if include_poles && !isempty(u.poles)
        return 0
    end
    prod_terms = []
    start = (u.special == true) ? 2 : 1
    for i = length(u):-1:start
        temp_sum = sum(u.freqs[length(u):-1:i])
        if !isequal(temp_sum, 0) push!(prod_terms, temp_sum) end
    end
    return isempty(prod_terms) ? 1 : prod(prod_terms)
end

Base.sum(u::BVector) = isempty(u) ? 0 : sum(u.freqs) 
Base.length(u::BVector) = length(u.freqs)
Base.getindex(u::BVector, i::Int) = u.freqs[i]
Base.size(u::BVector) = size(u.freqs)
