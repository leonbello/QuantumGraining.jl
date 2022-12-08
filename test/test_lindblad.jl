#Test file for lindblad.jl 
using QuantumCumulants

module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/lindblad.jl")
    include("../src/printing.jl")

    diagram = [(2,1), (1,0)]

    N = 2 #Number of atoms 
    h = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
    σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k)

    σz1 = 2*σ(1,1,1) - 1
    σz2 = 2*σ(1,1,2) - 1

    σx1 = σ(1,2,1) + σ(2,1,1)
    σx2 = σ(1,2,2) + σ(2,1,2)

    @cnumbers ω1 ω2 ωd g ϵd

    #Define ordered lists of frequencies and hamiltonian terms 
    ω_list = [0.001,0.001,0.001,ωd,-ωd]
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