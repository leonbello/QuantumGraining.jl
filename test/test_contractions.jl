using Revise
using QuantumCumulants
using QuantumGraining
using SymbolicUtils
using IterTools
using Symbolics
#using Test

#@testset "contractions" begin
module Tst
    using Test
    include("../src/diagrams.jl")
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
    μ1 = [0, 1, -1] # poles at 3
    ν1 = [5]
    ω = [(μ1, ν1)]
    diagram_correction(ω)


    ## OLD TESTS ##
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

    @cnumbers e1 e2 e3 f1 f2 f3
    ex = [(e1, f1),(e2,f2),(e3,f3)]
    l, r = ex
    @show l
    @show r

     @cnumbers a1 a2 a3 b1 b2 b3 b4 b5 b6 b7   
    diagram = [(1,2), (3,4)]
    ω = [a1,a2,a3,b1,b2,b3,b4,b5,b6,b7]
    μ, ν = [],[]
    ind = 0
    for (i, bubble) in enumerate(diagram)
        (μ_len, ν_len) = bubble
        push!(μ, ω[ind+1:ind+μ_len])
        ind += μ_len
        push!(ν, ω[ind+1:ind+ν_len])
        ind += ν_len
    end
end
@show ν

# if denominator goes to zero -> take limit


#=
end #testset
=#