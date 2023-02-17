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

    h_cav = FockSpace(:cavity)
    h_atom = NLevelSpace(:atom, (:g,:e))
    h = tensor(h_cav, h_atom)

    @qnumbers a::Destroy(h) σ::Transition(h)
    @cnumbers ω_0 ω_d g κ γ ϵ ω_3 ω_1

    σz = 2* σ(:e, :e)  - 1
    σx = σ(:g, :e) + σ(:e, :g)
    σy = 1im* (σ(:e, :g) - σ(:g, :e))
    @syms t::Real
    H_rab = ω_0 * a'*a + ω_d/2 * σz + g*(a' + a)*(σx) + ϵ*(a*exp(1im*ω_3*t) + a'*exp(-1im*ω_3*t))
    #@register f(t)
    ωs = [0,0,0,ω_3, -ω_3]
    hs = [ω_0 * a'*a,ω_d/2 * σz,g*(a' + a)*(σx),ϵ*a, ϵ*a']

    ωs2 = [ω_1,-ω_1,ω_3, -ω_3]
    hs2 = [0.5*ω_0 * a'*a,0.5*ω_0 * a'*a,ϵ*a, ϵ*a'] 
   

    k=1
    @cnumbers t
    H = simplify(effective_hamiltonian(k, ωs, hs))
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