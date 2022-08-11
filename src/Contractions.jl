using QuantumCumulants

"""
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

    dmodes = [Symbol(:ν, i) for i in range(1, m)]
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


