struct Diagram{T1, T2} <: AbstractVector{Tuple{T1, T2}}
    ω::Vector{Tuple{BVector{T1}, BVector{T2}}}
    up_poles::Vector{Tuple{Vector{Int}, Vector{Int}}}
    down_poles::Vector{Tuple{Vector{Int}, Vector{Int}}}

    function Diagram(ω::Vector{Tuple{BVector{T1}, BVector{T2}}}) where {T1,T2}
        up_poles, down_poles = find_all_poles(ω)
        return new{T1, T2}(ω, up_poles, down_poles)
    end
end
