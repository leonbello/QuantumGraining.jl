using Revise
using QuantumCumulants
using QuantumGraining

#using Test

#@testset "contractions" begin
module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")
    include("../src/expressions.jl")


    ## coeff ##
    c = bubble_coeff(3, 2)
    @definemodes μ 3
    @definemodes ν 2
    c(μ, ν)
    c([3, 2, 6], [1, 4, 7])
    typeof(c(μ,ν))

    c = bubble_coeff(5, 0)
    c(μ, ν)

    ## contraction_coeff ##
    c = contraction_coeff(2, 3)
    typeof(c(μ,ν))
    c(μ, ν)
    
    diagram = [(3, 0), (5, 0)]
    for (i, bubble) in enumerate(diagram)
        println("bubble $i: ($(bubble[1]), $(bubble[2])) ")
        println("----")
    end

    ## Hamiltonian contractions ##
    # check macro for global scope
    diagram2 = [(2,0), (3,0)]
    @definemodes x 3
    @show x
    @show μ

    @show get_max_modes(diagram)
    @show get_max_modes(diagram2)
    @definemodes x get_max_modes(diagram2)[1]
    @definemodes y get_max_modes(diagram2)[2]
    
    c1, cprod = calculate_coeff(x, y, τ, diagram2)
    @show c1
    @show cprod
    
    # global scope works fine
    m = 4
    @definemodes x m
    @show x

    typeof(diagram3)
    print(test_definemodes(diagram3))
end

# if denominator goes to zero -> take limit


#=
end #testset
=#