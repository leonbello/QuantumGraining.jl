#Test file for lindblad.jl 
using QuantumCumulants
using IterTools

module Tst
    using Test
    using IterTools
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/lindblad.jl")
    include("../src/printing.jl")
    include("../src/poles.jl")

    ## problem definition ##
    num_bubbles = 2
    @cnumbers a b c
    μ1 = [0, a, -a] # pole at 1 and 3
    ν1 = [a, -a, b, c] # pole at 2

    μ2 = [a, b, c, -(a+b+c)] # pole at 4
    ν2 = [3*c, 2*b] # no poles

    ω = [(μ1, ν1), (μ2, ν2)]

    s_list = [[1, 3], [4]]
    stag_list = [[2], []]                          # singular indices
    total_num_poles = 4
    
    sols = find_integer_solutions(3*num_bubbles, total_num_poles)     
    unl_list = reshape_sols(sols, total_num_poles, num_bubbles)           # partition for the inner sum

    ## singular_expansion & calculate_bubble_factor() - bubble 1 ##
    s = s_list[1]
    stag = stag_list[1]
    
    @show s, stag
    t1 = singular_expansion(μ1, ν1, unl_list[:, 1], s, stag)
    t2 = calculate_bubble_factor(ω, 1, unl_list[:, 1], total_num_poles, s, stag)
    
    ## diagram_correction(ω) ##
    t3 = diagram_correction(ω)

    # non-singular bubble
    μ1 = [1, 2, 3]     # 2*(2 + 3) -- 3*(2 + 3)
    #μ1 = [1, 3, 2]
    ν1 = [7, 4]
    ω = [(μ1, ν1)]
    diagram_correction(ω)
    
    
    μ1 = [0, 2, -2]
    ν1 = [4, -4]

    μ2 = [1, 3]
    ν2 = [5]
    ω = [(μ1, ν1), (μ2, ν2)]

    test = diagram_correction(ω)

    ## Old tests, may be irrelevant now ##
    diagram = [(2,1), (1,0)]

    N = 2 # number of atoms 
    h = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
    σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k)

    σz1 = 2*σ(1,1,1) - 1
    σz2 = 2*σ(1,1,2) - 1

    σx1 = σ(1,2,1) + σ(2,1,1)
    σx2 = σ(1,2,2) + σ(2,1,2)

    @cnumbers ω1 ω2 ωd g ϵd

    #Define ordered lists of frequencies and hamiltonian terms 
    ω_list = [0, 0, 0, ωd, -ωd]
    h_list = [ω1*σz1/2,ω2*σz2/2, g* σx1 * σx2, ϵd* σx2,ϵd* σx2]
    
    ## repeated_combinations ##
    ω_combos = repeated_combinations(ω_list, 5)
    unique(ω_combos)

    h_combos = repeated_combinations(h_list, 5)
    unique(simplify.(h_combos))                                 # seems to be an error with the simplification

    ## effective_hamiltonian ##
    k = 2
    ω_combos = repeated_combinations(ω_list, k)
    h_combos = repeated_combinations(h_list, k)
    h_eff = effective_hamiltonian(k, ω_list, h_list)

    # problem with 0 frequencies and limits
end