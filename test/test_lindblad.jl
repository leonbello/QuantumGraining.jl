#Test file for lindblad.jl 
using QuantumCumulants
using IterTools
using Symbolics

module Tst
    using Test
    using IterTools
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/lindblad.jl")
    include("../src/printing.jl")
    include("../src/poles.jl")



    """
    ## Old tests, may be irrelevant now ##
    """
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
    
    
    ω2 = [([1, ωd], [1,1, -ωd])]
    @show diagram_correction(ω2)

    ## repeated_combinations ##
    ω_combos = repeated_combinations(ω_list, 5)
    unique(ω_combos)

    h_combos = repeated_combinations(h_list, 5)
    unique(simplify.(h_combos))                                 # seems to be an error with the simplification

    ## effective_hamiltonian ##
    k = 1
    ω_combos = repeated_combinations(ω_list, k)
    h_combos = repeated_combinations(h_list, k)
    h_eff = effective_hamiltonian(k, ω_list, h_list)
    #simplify(h_eff)
    

    # problem with 0 frequencies and limits
end