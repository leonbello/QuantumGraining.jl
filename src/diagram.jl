"""
    Diagram{T1, T2} <: AbstractVector{Tuple{T1, T2}}
A struct representing a diagram, also takes care of the edge mode.

Arguments:
    - bubbles: A vector of bubbles.
    - freqs: A vector of tuples of BVectors. The first element of the tuple is the up-bubble, the second is the down-bubble.
    - shape: A vector of tuples representing the shape of the bubbles.
    - up_poles: A vector of vectors representing the up-poles of the bubbles.
    - down_poles: A vector of vectors representing the down-poles of the bubbles.
    - num_poles: An integer representing the number of poles in the diagram.
"""
struct Diagram{T1, T2} <: AbstractVector{Tuple{T1, T2}}
    bubbles::Vector{Bubble{T1, T2}}
    freqs::Vector{Tuple{BVector{T1}, BVector{T2}}}
    shape::Vector{Tuple{Int, Int}}
    up_poles::Vector{Vector{Int}}
    down_poles::Vector{Vector{Int}}
    num_poles::Int

    """
        Diagram(freqs::Vector{Tuple{BVector{T1}, BVector{T2}}}
    A Diagram constructor taking a vector of tuples of BVectors. Useful if you want to specifically define the special mode manually.

    Arguments:
        - freqs: A vector of tuples of BVectors. The first element of the tuple is the up-bubble, the second is the down-bubble.
    Returns:
        - Diagram: A Diagram struct. Note that the first bubble must be denoted special for normal operation.
    """
    function Diagram(bubbles::Vector{Bubble{T1, T2}}) where {T1, T2}
        freqs = [(b.up, b.down) for b in bubbles]
        up_poles, down_poles = find_all_poles(freqs)
        shape = [(length(μ), length(ν)) for (μ, ν) in freqs]
        num_poles = count_poles(up_poles, down_poles)
        return new{T1, T2}(bubbles, freqs, shape, up_poles, down_poles, num_poles)
    end
    
    """
        Diagram(bubbles::Vector{Bubble{T1, T2}})
    A Diagram constructor taking a vector of bubbles. Useful if you want to specifically define the special bubble manually.
    Arguments:
        - bubbles: A vector of bubbles.
    Returns:
        - Diagram: A Diagram struct. Note that the first bubble must be denoted special for normal operation.
    """
    function Diagram(freqs::Vector{Tuple{BVector{T1}, BVector{T2}}}) where {T1,T2}
        bubbles = Vector{Bubble{T1, T2}}()
        for (μ, ν) in freqs
            special = μ.special || ν.special
            push!(bubbles, Bubble(μ, ν, special))
        end
        return Diagram(bubbles)
    end 

    """
        Diagram(freqs::Vector{Tuple{Vector{T1}, Vector{T2}}})
    A Diagram constructor taking a vector of tuples of vectors. Useful if you want to automatically define the special bubble.

    Arguments:
        - freqs: A vector of tuples of vectors. The first element of the tuple hold the up-bubble frequencies, the second is the down-bubble frequencies.
    Returns:
        - Diagram: A Diagram struct. Note that the first bubble is denoted a special edge mode bubble automatically.
    """
    function Diagram(freqs::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1,T2}
        bubbles = Vector{Bubble{T1, T2}}()
        for (idx, (μ, ν)) in enumerate(freqs)
            push!(bubbles, Bubble(UVec(μ), DVec(ν), special=(idx == 1)))
        end
        return Diagram(bubbles)
    end
end
Base.size(d::Diagram) = size(d.freqs)
Base.getindex(d::Diagram, i::Int) = d.bubbles[i]
Base.length(d::Diagram) = length(d.freqs)
Base.push!(d::Diagram, freq::Tuple{BVector{T1}, BVector{T2}}) where {T1,T2} = push!(d.freqs, freq)

