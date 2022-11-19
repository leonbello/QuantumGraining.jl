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
    _definemodes(name, n)

Returns an expression of the form @cnumbers name1 name2...
"""
function _definemodes(name::Base.Symbol, n::Int, fmt=:QuantumCumulants)
    new_ex = :(@cnumbers)
    sym_arr = []
    for i in 1:n
        ex = Symbol(name, i)
        sym_arr = push!(sym_arr, ex)
        push!(new_ex.args, ex)
    end
    return (new_ex, sym_arr)
end

"""
    @definemodes

Helper macro to define `n` up-modes and `m` down-modes. Generates two vectors of symbols.
        
# Examples
    @definemodes μ 3 ν 5    
    > μ = [μ1 μ2 μ3]
    > ν = [ν1 ν2 ν3 ν4 ν5]
    > typeof(μ) = Array{QuantumCumulants.cnumber}
"""
macro definemodes(uname, unum)
    uexp, u_arr = _definemodes(uname, eval(unum))
    esc_uname = esc(uname)
    return quote
        @cnumbers τ
        $uexp
        $esc_uname = eval.($(u_arr))
    end
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

"""
    coeff(μ, ν, endmode)

    Calculates the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient.

    Arguments:
    - μ: a vector of the symbolic up-mode frequencies in the bubble.
    - ν: a vector of the symbolic down-mode frequencies in the bubble.

    Returns:
    - C_{l,r}(μ, ν): a symbolic expression for the bubble coeffeicient.
"""
function coeff(μ, ν, τ, endmode=0) 
    l = length(μ)
    r = length(ν)
    sum_μ = l > 0 ? sum(μ) : 0
    sum_ν = r > 0 ? sum(ν) : 0
    f(ω, τ) = exp(-1/2*ω^2*τ^2)                      # define filter function
    s = (-1)^(2*l + r)                               # symmetry sign factor
    #in general (-1)^(num_bubbles + l)
    return s*f(sum_μ + sum_ν, τ)/(_Ω(μ[1:(end-endmode)])*_Ω(ν)) 
end

"""
    Given a diagram, returns a vector of the terms with the coeffcients for the contraction expression.
"""
function calculate_coeff(μ, ν, τ, diagram::Array{Tuple{Int64, Int64}})
    h = []
    endmode = 0
    for (i, bubble) in enumerate(diagram)
        if i == length(diagram)
            endmode = 1  
        end
        n, m = bubble
        push!(h, coeff(μ[1:n], ν[1:m], τ, endmode))
    end
    return prod(h), h
end



# """
#     Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
# """
# function coeff(n, m, endmode=0) 
#     f(ω, τ) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
#     return f(sum(μ) + sum(ν), τ)/(_Ω(μ)*_Ω(ν)) 
# end

#####
# """
#     Contraction

# A struct representing a contraction, includes all terms.

# # Fields
# - `order::Tuple{Int, Int}`: A tuple representing the order of the contraction.
# - `root::AbstractDiagramNode`: Root node of the decomposition tree.
# - `diagrams::Array`: An array holding all the diagrams. 
# """
# struct Contraction
#     order::Tuple{Int, Int}
#     root::DiagramNode
#     diagrams::Array

#     function Contraction(order)
#         order = order
#         root = DiagramNode(order)
#         diagrams = get_diagrams(root)
#     end
# end

# macro definemodes(n, m) 

#     umodes = eval(:([Symbol(:μ, i) for i in range(1, $n)]))
#     ex1 = :(@cnumbers)
#     for mode in umodes
#         push!(ex1.args, mode)
#     end                                                                     # @cnumbers μ1 μ2 ... μn

#     dmodes = [Symbol(:ν, i) for i in range(1, eval(:($m)) )]
#     ex2 = :(@cnumbers)
#     for mode in dmodes
#         push!(ex2.args, mode)
#     end                                                                     # @cnumbers ν1 ν2 ... νm
    
#     ex3 = quote
#         μ = [eval(Symbol(:μ, i)) for i in range(1, $n)]
#         ν = [eval(Symbol(:ν, i)) for i in range(1, $m)]
#     end                                                                     # generates a list of the mode variables

#     ex4 = quote
#         @cnumbers τ  
#     end
#     return esc(:($ex1; $ex2; $ex3; $ex4))
# end

# macro definemodes(sym1, sym2, diagram)
#     #sym1 = esc(sym1)
#     #sym2 = esc(sym2)                                                        # take the names the user gives
#     n, m = eval(:(_maxmodes($diagram)))

#     umodes = [Symbol(sym1, i) for i in range(1, n)]
#     #umodes = eval( :([Symbol(sym1, i) for i in range(1, $n)]) )
#     ex1 = :(@cnumbers)
#     for mode in umodes
#         push!(ex1.args, mode)
#     end                                                                     # @cnumbers μ1 μ2 ... μn

#     dmodes = [Symbol(sym2, i) for i in range(1, m)]
#     ex2 = :(@cnumbers)
#     for mode in dmodes
#         push!(ex2.args, mode)
#     end                                                                     # @cnumbers ν1 ν2 ... νm
    
#     μ = Symbol(sym1, "_arr")
#     ν = Symbol(sym2, "_arr")
#     ex3 = quote
#         $(μ) = [eval(Symbol($sym1, i)) for i in range(1, $n)]
#         $(ν) = [eval(Symbol($sym2, i)) for i in range(1, $m)]
#     end                                                                     # generates a list of the mode variables

#     ex4 = quote
#         @cnumbers τ  
#     end
#     return esc(:($ex1; $ex2; $ex3; $ex4))
# end