using QuantumCumulants
using QuantumGraining

#using Test

#@testset "contractions" begin
module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")

    diagram = [(3, 4), (5, 1)]
    for (i, bubble) in enumerate(diagram)
        println("bubble $i: ($(bubble[1]), $(bubble[2])) ")
        println("----")
    end

    n, m = get_max_modes(diagram)
    @definemodes μ n
    @definemodes ν m
    c = coeff(μ, ν, τ)
    C, h = calculate_coeff(μ, ν, τ, diagram)
    @show h[1]
    @show h[2]
    @show C
end

# if denominator goes to zero -> take limit


#=
end #testset
=#