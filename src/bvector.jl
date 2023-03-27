struct BVector{T} <: AbstractVector{T}
    vec::Vector{T}
    poles::Vector{Int}
    special::Bool
    type::Symbol

    function BVector(vec::Vector{T}, type::Symbol, special::Bool) where {T}
        if type == :up
            vec = reverse(vec)
        elseif type == :down
            vec = vec
        else
            error("type must be :up or :down")
        end
        poles = find_poles(vec)
        return new{T}(vec, poles, special, type)
    end
end

function UVec(u::Vector; special=false)
    return BVector(u, :up, special)
end

function DVec(u::Vector; special=false)
    return BVector(u, :down, special)
end

"""
    vec_factorial(u; include_poles=true)

Calculates vector factorial for a frequency list u. 
"""
function vec_factorial(u; special = false, include_poles=true)
    
    if isempty(u)
        return 1
    end

    if special deleteat!(u, 1) end

    prod_terms = []
    for i = 1:length(u)
        temp_sum = sum(u[1:i])
        if isequal(temp_sum, 0)
            if include_poles == true
                temp_prod = 0
            end
        else
            push!(prod_terms, temp_sum)
        end
    end
    return isempty(prod_terms) ? 1 : prod(prod_terms)
end

find_poles(u::BVector) = find_poles(u.vec)
function find_all_poles(ω::Vector{Tuple{BVector{T1}, BVector{T2}}}) where {T1, T2}
    ω = [(μ.vec, ν.vec) for (μ, ν) in ω]
    return find_all_poles(ω)
end 

Base.sum(u::BVector) = isempty(u) ? 0 : sum(u.vec) 
Base.factorial(u::BVector) = vec_factorial(u.vec, u.special, include_polse=false) # TODO: include poles
Base.length(u::BVector) = length(u.vec)
Base.getindex(u::BVector, i::Int) = u.vec[i]
Base.size(u::BVector) = size(u.vec)

struct Diagram{T1, T2} <: AbstractVector{Tuple{T1, T2}}
    freqs::Vector{Tuple{BVector{T1}, BVector{T2}}}
    shape::Vector{Tuple{Int, Int}}
    up_poles::Vector{Vector{Int}}
    down_poles::Vector{Vector{Int}}
    num_poles::Int
    
    function Diagram(freqs::Vector{Tuple{BVector{T1}, BVector{T2}}}) where {T1,T2}
        up_poles, down_poles = find_all_poles(freqs)
        shape = [(length(μ), length(ν)) for (μ, ν) in freqs]
        num_poles = count_poles(up_poles, down_poles)
        return new{T1, T2}(freqs, shape, up_poles, down_poles, num_poles)
    end
end
Base.size(d::Diagram) = size(d.freqs)
Base.getindex(d::Diagram, i::Int) = d.freqs[i]
Base.length(d::Diagram) = length(d.freqs)