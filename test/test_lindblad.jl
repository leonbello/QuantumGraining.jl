#Test file for lindblad.jl 
using QuantumCumulants
using IterTools
using Symbolics
using SymbolicUtils.Rewriters

module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/lindblad.jl")
    include("../src/printing.jl")
    include("../src/poles.jl")

    #=
        Old tests, may be irrelevant now ##
    =#
    diagram = [(2,1), (1,0)]

    
    N = 2 # number of atoms 
    h = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
    σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k)

    σz1 = 2*σ(1,1,1) - 1
    σz2 = 2*σ(1,1,2) - 1

    σx1 = σ(1,2,1) + σ(2,1,1)
    σx2 = σ(1,2,2) + σ(2,1,2)

    #@show typeof(σz1)
    @cnumbers ω_1 ω_2 ω_d g ϵ_d
    typeof(ω_1)
    #Define ordered lists of frequencies and hamiltonian terms 
    ωs = [0, 0, 0, ω_d, -ω_d]
    hs = [ω_1*σz1/2,ω_2*σz2/2, g* σx1 * σx2, ϵ_d* σx2,ϵ_d* σx2]
    typeof(hs[2])

    k=1
    @cnumbers t
    H = simplify(effective_hamiltonian(k, ωs, hs,t))
    #collect_rule = @rule(+(~~xs) => ~~xs)
    #H_list = simplify(collect_rule(H))
    #typeof(ω_d*σz1)
    secondOrderRule = @acrule((~a)*(~y)*(~z) + (~b)*(~y)*(~z) => ((~a)+(~b))*(~y)*(~z))
    firstOrderRule = @acrule((~a)*(~y) + (~b)*(~y) => ((~a)+(~b))*(~y))
    rules = Chain([firstOrderRule,secondOrderRule])
    H1 = Fixpoint(rules)(H)
    #r = Chain[]

    typeof(H1)
    eqs = meanfield([σz1], H1, []; order=1)
    eqs1 = complete(eqs)
    ω2 = [([1, ω_d], [1,1, -ω_d])]
    @show diagram_correction(ω2)
    ω3 = [1,1,1,ω_d, -ω_d]
    @show contraction_coeff((3,0), ω3) + contraction_coeff((0,3), -reverse(ω3))


    ## repeated_combinations ##
    #ω_combos = repeated_combinations(ω_list, 5)
    #unique(ω_combos)

    #h_combos = repeated_combinations(h_list, 5)
    #unique(simplify.(h_combos))                                 # seems to be an error with the simplification

    ## effective_hamiltonian ##
    #k = 2
    #ω_combos = repeated_combinations(ω_list, k)
    #h_combos = repeated_combinations(h_list, k)
    #@cnumbers t
    #h_eff = effective_hamiltonian(k, ω_list, h_list,t)
    #simplify(h_eff)
    


    # problem with 0 frequencies and limits
end