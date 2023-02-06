using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
#using Test

#@testset "contractions" begin
module Tst
    using Test
    using IterTools
    include("../src/diagrams.jl")
    include("../src/poles.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")
    include("../src/expressions.jl")

    ## problem definition ##
    num_bubbles = 2
    @cnumbers a b c
    μ1 = [0, a, -a] # pole at 1 and 3
    ν1 = [a, -a, b, c] # pole at 2

    μ2 = [a, b, c, -(a+b+c)] # pole at 4
    ν2 = [3*c, 2*b] # no poles

    ω = [(μ1, ν1), (μ2, ν2)]

    s_list = [[1, 3], [4]]
    stag_list = [[2], []]  # singular indices
    total_num_poles = 4
    
    num_vars = 3*num_bubbles
    num_sols = binomial(total_num_poles + num_vars - 1, num_vars - 1)
    sols = find_integer_solutions(3*num_bubbles, total_num_poles)     
    unl_list = reshape_sols(sols, total_num_poles, num_bubbles)           # partition for the inner sum

    ## singular_expansion & calculate_bubble_factor() - bubble 1 ##
    s = s_list[1]
    stag = stag_list[1]
    
    @show s, stag
    t1 = singular_expansion(μ1, ν1, unl_list[:, 1], s, stag)
    t2 = calculate_bubble_factor(ω, 1, unl_list[:, 1], s, stag)
    
    ## diagram_correction(ω) ##
    t3 = diagram_correction(ω)

    # non-singular bubble
    μ1 = [1, 2, 3]     # 2*(2 + 3) -- 3*(2 + 3)
    #μ1 = [1, 3, 2]
    ν1 = [7, 4]
    ω = [(μ1, ν1)]
    diagram_correction(ω)
    
    μ1 = [1, 2, 3]
    ν1 = [7, 4]

    μ2 = [1, 3]
    ν2 = [5]
    ω = [(μ1, ν1), (μ2, ν2)]
    test = diagram_correction(ω)

    # only up-bubbles
    μ1 = [1, 4]
    ν1 = []
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    μ1 = []
    ν1 = [1, 4]
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    # singular bubbles
    μ1 = [0, 1, -1] # pole at 2
    ν1 = [5]
    ω = [(μ1, ν1)]

    diagram_correction(ω)
end


#=
end #testset
=#