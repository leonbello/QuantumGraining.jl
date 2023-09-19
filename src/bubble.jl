using QuantumGraining

"""
    Bubble{T1, T2} <: AbstractVector{Tuple{T1, T2}}    
A struct representing a bubble, also takes care of the edge case.

Arguments:
    - up: A BVector representing the up-bubble.
    - down: A BVector representing the down-bubble.
    - shape: A tuple representing the shape of the bubble.
    - special: A boolean representing whether the bubble is a special edge mode bubble.
"""
mutable struct Bubble{T1, T2} <: AbstractVector{Tuple{T1, T2}}
    up::BVector{T1}
    down::BVector{T2}
    shape::Tuple{Int, Int}
    freqs::Vector{Union{T1, T2}}
    special::Bool

    """
        Bubble(freqs::Vector{T}, shape::Tuple{Int, Int}, special::Bool)
    A Bubble constructor taking a shape `(n, m)` and a vector of frequencies in the form [μ1, μ2, μn, ν1, ν2, νm]. 

    Arguments:
        - freqs: A vector of tuples of vectors. The first element of the tuple is the up-bubble, the second is the down-bubble.
        - shape: A tuple representing the shape of the bubble.
        - special: A boolean representing whether the bubble is a special edge mode bubble.
    """
    function Bubble(freqs::Vector{T}, shape::Tuple{Int, Int}; special = false) where {T}
        if length(freqs) != (shape[1] + shape[2])
            error("Shape does not match vector length")
        end
        up = UVec(freqs[1:shape[1]]; special=special)
        down = DVec(freqs[(shape[1] + 1):(shape[1] + shape[2])]; special=special)
        return new{T, T}(up, down, shape, freqs, special)
    end

    """
        Bubble(up::BVector{T1}, down::BVector{T2}, special::Bool)
    A Bubble constructor taking two BVectors. Useful if you want to automatically define the special bubble manually.

    Arguments:
        - up: A BVector representing the up-bubble.
        - down: A BVector representing the down-bubble.
        - special: A boolean representing whether the bubble is a special edge mode bubble.
    """
    function Bubble(up::BVector{T1}, down::BVector{T2}; special::Bool=false) where {T1,T2}
        if isempty(up) && isempty(down)
            error("Bubble must have at least one frequency")
        end

        if special
            if isempty(up) && isempty(down)
                error("Special bubble must have at least one frequency")
            elseif isempty(up) && !isempty(down)
                up.special = false
                down.special = true
            elseif !isempty(up)
                up.special = true
                down.special = false
            end
        end
        
        freqs = [up.freqs..., reverse(down.freqs)...]
        shape = (length(up), length(down))
        return new{T1, T2}(up, down, shape, freqs, special)
    end


    """
        Bubble(up::Vector{T1}, down::Vector{T2}, special::Bool)
    A Bubble constructor taking two vectors. Useful if you want to automatically define the special bubble manually.
    Arguments:
        - up: A vector representing the up-bubble.
        - down: A vector representing the down-bubble.
        - special: A boolean representing whether the bubble is a special edge mode bubble.
    """
    function Bubble(up::Vector{T1}, down::Vector{T2}; special::Bool=false) where {T1, T2}
        return Bubble(UVec(up,special=special), DVec(down); special=special)
    end
end
Base.length(b::Bubble) = prod(b.shape)

function Base.getindex(b::Bubble, i::Int)
    vec = [b.up..., reverse(b.down)...]
    return vec[i]
end
Base.size(b::Bubble) = b.shape