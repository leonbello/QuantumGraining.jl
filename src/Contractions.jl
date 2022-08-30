using QuantumCumulants

"""
<<<<<<< HEAD
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
    =#
    #=
    Macros pass the symbols into the macro code, not the variables. 
    Anywhere we write `n` or `m` it would literally be the symbols, we usually avoid that by writing $(var"name") inside
    the expressions, but it cannot be used outside of the expression context.
    Since some of the expression are evaluated directly, we use eval(var"name").
    =#
    umodes = [Symbol(:μ, i) for i in range(1, eval(:($(n))) )]
    ex1 = :(@cnumbers)
    for mode in umodes
        push!(ex1.args, mode)
    end                                                                     # @cnumbers μ1 μ2 ... μn

    dmodes = [Symbol(:ν, i) for i in range(1, eval(:($m)) )]
=======
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
>>>>>>> contraction-coeffecients
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

<<<<<<< HEAD

"""
    Simple helper function to calculate the normalization for the correction terms.
"""
function _Ω(freqs)
    n = length(freqs)
    terms = [sum(freqs[i:end]) for i in 1:n]
    return prod(terms)
end

"""
    Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
"""
function coeff(μ, ν, endmode=0) 
    f(ω) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
    return f(sum(μ) + sum(ν))/(_Ω(μ)*_Ω(ν)) 
end

function coeff(bubble, endmode=0)
    coeff(μ[1:bubble[1]], ν[1:bubble[2]], endmode)
end
=======
"""
    Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
"""
function coeff(n, m, endmode=0) 
    f(ω, τ) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
    return f(sum(μ) + sum(ν), τ)/(_Ω(μ)*_Ω(ν)) 
end

>>>>>>> contraction-coeffecients

"""
    Given a diagram, returns a vector of the terms with the coeffcients for the effective Hamiltonian expression.
"""
function calculate_coeff(diagram::Array{Tuple{Int64, Int64}})
<<<<<<< HEAD
    h = []
    ububs = zeros(Int64, length(diagram))
    dbubs = zeros(Int64, length(diagram))
    for (i, bubble) in enumerate(diagram)
        ububs[i] = bubble[1]
        dbubs[i] = bubble[2]
    end
    umax = maximum(ububs)
    dmax = maximum(dbubs)                            # to avoid redundancy, we find the maximum amount of participating modes
    @show umax
    @show dmax
    #@definemodes umax dmax                       # define all the participating modes as symbols

    for (i, bubble) in enumerate(diagram)
        push!(h, coeff(bubble))
=======

    h = []
    for (i, bubble) in enumerate(diagram)
        @definemodes bubble[1] bubble[2]
        push!(h, _calculate_coeff(bubble))
>>>>>>> contraction-coeffecients
    end
    return h
end


<<<<<<< HEAD
# Testing area
diagram = [(3, 4), (2, 1)]
for (i, bubble) in enumerate(diagram)
    @show bubble[1]
    @show bubble[2]
    println("----")
end

@definemodes diagram[1][1] diagram[1][2]
@show μ
@show ν
@show τ

@definemodes diagram[2][1] diagram[2][2]
@show μ
@show ν
@show τ

umax = 3
dmax = 4
@definemodes umax dmax
@show μ
@show ν
@show τ

c = calculate_coeff(diagram)
@show c


=======
>>>>>>> contraction-coeffecients
