using QuantumCumulants

"""
    Contraction

A struct representing a contraction.

# Fields
- `order::Tuple{Int, Int}`: A tuple representing the order of the contraction.
- `root::AbstractDiagramNode`: Root node of the decomposition tree.
- `diagrams::Array`: An array holding all the diagrams. 
"""
struct Contraction
    order::Tuple{Int, Int}
    root::DiagramNode
    diagrams::Array

    function Contraction(order)
        order = order
        root = DiagramNode(order)
        diagrams = get_diagrams(root)
    end
end

"""
    @definemodes

Helper macro to define `n` up-modes and `m` down-modes. Generates two vectors of symbols.
        
# Examples
    @definemodes 3 5    
    > μ = [μ1 μ2 μ3]
    > ν = [ν1 ν2 ν3 ν4 ν5]
    > typeof(μ) = QuantumCumulants.cnumber
"""
macro definemodes(n, m) 
    #=  
        TODO:
        - Change the macro so the user can define their own variable names
        - Test; there are some weird undefined variable errors
    =#
    umodes = eval(:([Symbol(:μ, i) for i in range(1, $n)]))
    ex1 = :(@cnumbers)
    for mode in umodes
        push!(ex1.args, mode)
    end                                                                     # @cnumbers μ1 μ2 ... μn

    dmodes = [Symbol(:ν, i) for i in range(1, eval(:($m)) )]
"""
    _Ω(freqs...)
Simple helper function to calculate the normalization for the correction terms.
"""
function _Ω(freqs)
    n = length(freqs)
    terms = [sum(freqs[i:end]) for i in 1:n]
    return prod(terms)
end

"""
    Helper macro to define `n` up-modes and `m` down-modes
"""
macro definemodes(n, m)
    umodes = [Symbol(:μ, i) for i in range(1, n)]
    ex1 = :(@cnumbers)
    for mode in umodes
        push!(ex1.args, mode)
    end                                                                      # @cnumbers μ1 μ2 ... μn

    dmodes = eval(:([Symbol(:ν, i) for i in range(1, m)]))
    ex2 = :(@cnumbers)
    for mode in dmodes
        push!(ex2.args, mode)
    end                                                                     # @cnumbers ν1 ν2 ... νm
    
    ex3 = quote
        μ = [eval(Symbol(:μ, i)) for i in range(1, $n)]
        ν = [eval(Symbol(:ν, i)) for i in range(1, $m)]
    end                                                                     # generates a list of the mode variables

    ex4 = quote
        @cnumbers τ  
    end
    return esc(:($ex1; $ex2; $ex3; $ex4))
end


"""
    Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
"""
function coeff(μ, ν, endmode=0) 
    f(ω) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
    return f(sum(μ) + sum(ν))/(_Ω(μ[1:(end-endmode)])*_Ω(ν)) 
end

#=
function coeff(bubble, endmode=0)
    # add symmetry sign factors
    coeff(μ[1:bubble[1]], ν[1:bubble[2]], endmode)
end
=#

"""
    Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
"""
function coeff(n, m, endmode=0) 
    f(ω, τ) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
    return f(sum(μ) + sum(ν), τ)/(_Ω(μ)*_Ω(ν)) 
end

"""
    Given a diagram, returns a vector of the terms with the coeffcients for the effective Hamiltonian expression.
"""
function calculate_coeff(diagram::Array{Tuple{Int64, Int64}})
    h = []
    for (i, bubble) in enumerate(diagram)
        @definemodes bubble[1] bubble[2]
        push!(h, _calculate_coeff(bubble))
    end
    return h
end

function _maxmodes(diagram::Array{Tuple{Int64, Int64}})
    ububs = zeros(Int64, length(diagram))
    dbubs = zeros(Int64, length(diagram))
    for (i, bubble) in enumerate(diagram)
        ububs[i] = bubble[1]
        dbubs[i] = bubble[2]
    end
    return (maximum(ububs), maximum(dbubs))
end
