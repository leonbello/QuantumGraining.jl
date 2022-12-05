using QuantumCumulants
using SymbolicUtils
"""
    _Ω(freqs...)
Simple helper function to calculate the normalization factor for the correction terms.
"""
function _Ω(freqs)
    n = length(freqs)
    terms = [sum(freqs[i:end]) for i in 1:n]
    return n > 0 ? prod(terms) : 1
end

"""
    coeff(l, r, endmode, freqs...)

    Calculates the coefficient of a single bubble, 
    returns a symbolic function that evaluates to the coefficient as a function of frequency.

    Arguments:
    - μ: a vector of the symbolic up-mode frequencies in the bubble.
    - ν: a vector of the symbolic down-mode frequencies in the bubble.

    Returns:
    - C_{l,r}(μ, ν): a symbolic expression for the bubble coeffeicient.
"""
function coeff(l, r, endmode=0) 
    @cnumbers τ
    f(ω, τ) = exp(-1/2*ω^2*τ^2)                      # define filter function
    s = (-1)^(2*l + r)                               # symmetry sign factor, in general (-1)^(num_bubbles + l)

    # define an expression for the bubble with the frequencies (μ, ν) as function arguments
    function _c(μ, ν) 
        sum_μ = l > 0 ? sum(μ) : 0
        sum_ν = r > 0 ? sum(ν) : 0
        return s*f(sum_μ + sum_ν, τ)/(_Ω(μ[1:(end-endmode)])*_Ω(ν))
    end
    return _c                                                                    # should return a symbolic *function*
end

function contraction_coeff(left::Int, right::Int)
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)
    function _c(μ, ν)
        ex = 0                                                              # initialize an empty expression
        for diagram in diagrams
            diag_ex = 1
            for (i, bubble) in enumerate(diagram)
                endmode = (i == length(diagram)) ? 1 : 0
                diag_ex = diag_ex*(coeff(left, right, endmode))(μ, ν)        # calculate the bubble correction term and multiply it in
            end
            ex = ex + diag_ex                                   # add the diagram to the contraction coefficient expression
        end
        return ex
    end
    return _c
end

"""
    get_max_modes(diagram)

Returns the maximal number of modes in a diagram.
        
# Examples
    diagram = [(3,4), (2, 5)]
    n, m = get_max_modes(diagram)  
    > n = 3
    > m = 5
"""
function get_max_modes(diagram::Array{Tuple{Int64, Int64}})
    ububs = zeros(Int64, length(diagram))
    dbubs = zeros(Int64, length(diagram))
    for (i, bubble) in enumerate(diagram)
        ububs[i] = bubble[1]
        dbubs[i] = bubble[2]
    end
    return (maximum(ububs), maximum(dbubs))
end

# """
#     Given a diagram, returns a vector of the terms with the coeffcients for the contraction expression.
# """
# function calculate_coeff(μ, ν, τ, diagram::Array{Tuple{Int64, Int64}})
#     h = []
#     endmode = 0
#     for (i, bubble) in enumerate(diagram)
#         if i == length(diagram)
#             endmode = 1  
#         end
#         n, m = bubble
#         push!(h, coeff(μ[1:n], ν[1:m], τ, endmode))
#     end
#     return prod(h), h
# end

# """
#     Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
# """
# function coeff(n, m, endmode=0) 
#     f(ω, τ) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
#     return f(sum(μ) + sum(ν), τ)/(_Ω(μ)*_Ω(ν)) 
# end
