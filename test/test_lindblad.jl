#Test file for lindblad.jl 
using QuantumCumulants
using IterTools
using Symbolics
using SymbolicUtils.Rewriters
using SymPy
using SymbolicUtils
using Latexify

module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/bvector.jl")
    include("../src/bubble.jl")
    include("../src/diagram.jl")
    include("../src/poles.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")
    include("../src/corrections.jl")
    include("../src/lindblad.jl")

   
    
    #Setup
    h_cav = FockSpace(:cavity)
    h_atom = NLevelSpace(:atom, (:g,:e))
    h = tensor(h_cav, h_atom)

    @qnumbers a::Destroy(h) σ::Transition(h)
    SymPy.@syms ω_0 ω_d g κ γ ϵ ω_3 ω_1 
    SymPy.@syms ω_c::real ω_a::real

    σz = 2* σ(:e, :e)  - 1
    σx = σ(:g, :e) + σ(:e, :g)
    σy = 1im* (σ(:e, :g) - σ(:g, :e))
    SymPy.@syms t::Real


    #Symbolic test: Rabi Hamiltonian
    
    ωs_rab = [ω_c + ω_a, -ω_c - ω_a, -ω_c + ω_a, ω_c - ω_a]
    hs_rab = [a*σ(:e, :g),a'*σ(:g, :e),a*σ(:g, :e), a'*σ(:e, :g)]
    
    
    eff_ham_rab, ops_eff_rab, ωs_eff_rab, gs_eff_rab = effective_hamiltonian(hs_rab, ωs_rab, 2)
    sum(eff_ham_rab)
    ops_eff_rab
    render(latexify(to_symbol(gs_eff_rab[14])))
    
    #Numerical test
    ω_n = 1
    ω_m = 1

    
    ωs_rab_num = [ω_n + ω_m, -ω_n - ω_m, -ω_n + ω_m, ω_n - ω_m]
    hs_rab_num = [a*σ(:e, :g),a'*σ(:g, :e),a*σ(:g, :e), a'*σ(:e, :g)]
    
    eff_ham_rab_num, ops_eff_rab_num, ωs_eff_rab_num, gs_eff_rab_num = effective_hamiltonian(hs_rab_num, ωs_rab_num, 2)

    gs_eff_rab_num[2].polys
    render(latexify(to_symbol(gs_eff_rab_num[2])))
    eff_ham_rab_num

    to_symbol(contraction_coeff(2,0,[ω_a,ω_c]))
    to_symbol(contraction_coeff(0,2,[ω_a,ω_c]))
    """
    H_rab = ω_0 * a'*a + ω_d/2 * σz + g*(a' + a)*(σx) + ϵ*(a*exp(1im*ω_3*t) + a'*exp(-1im*ω_3*t))
    #@register f(t)
    a*a'
    ωs = [1,1,2,2.0, -2.0]
    hs = [ω_0 * a'*a,ω_d/2 * σz,0.5,ϵ*a, ϵ*a']
    #hs = [1,1,1,1,1]
    ops_eff, ωs_eff, gs_eff =  effective_hamiltonian(hs, ωs, 2)
    gs_eff
    eff_ham = 0
    for i in length(gs_eff)
        prefacs = gs_eff[i].prefacs
        expons = gs_eff[i].exponents
        eff_ham += prefacs[1]*exp(1im*expon*t)*(hs_eff[i]*exp(1im *sum(ωs_eff[i])*t))
    end

    perm_h, perm_Ω = repeated_combinations(hs, ωs,2)
    l = perm_h[1]
    typeof(l)
    #contraction_coeff(3, 0,l)
    contraction_coeff(2, 0,[1,1]) + contraction_coeff(2, 0,[2,1])
    effective_hamiltonian(hs,ωs,2)
    ωs2 = [ω_1,-ω_1,ω_3, -ω_3]
    hs2 = [0.5*ω_0 * a'*a,0.5*ω_0 * a'*a,ϵ*a, ϵ*a'] 
    example = effective_hamiltonian(hs2,ωs2,2)[1]
    sum(hs2)
    hs2_sym = Sym(sum(hs2))
    Sym(a'*a)
    collect(hs2_sym)
    
    eff_ham = effective_hamiltonian(hs2,ωs2,2)
    SymPy.@syms t::Real
    exp(complex(0,1))
    exp(-1im*sum(eff_ham[2][1]))
    
    using Latexify
    sum1 = 0
    for i in 1:4
        sum1 += to_symbol(eff_ham[3][1])*eff_ham[1][1]
    end
    @variables t::Real
    SymPy.@syms x y
    z = x + 1im*y
    exp(z)
    #Extract actual effective Hamiltonian
    eff_haml = 0
    for i in 1:length(eff_ham[1])
        eff_haml += exp(-1im*sum(eff_ham[2][i]))*to_symbol(eff_ham[3][i])*eff_ham[1][i]   
    end
    
    SymPy.simplify(eff_haml)
    arguments(eff_haml)
    operation(expr)
    eff_ham[3]
    
    # Define a function to split an expression by '+'
    split_expression(expr) = split(expr, '+')

    # Split each expression and flatten the resulting list
    split_elements = vcat(map(split_expression, example)...)


    # Define symbolic variables
    @syms x y z

    # Create a symbolic expression
    expr = x + 2*y - 3*z^2

    # Extract all the components from the expression
    components = args(expr)





    #Contraction coefficient test
    ωs = [1,1,0,1, -2.0]
    contraction_coeff(2,0,ωs2) + contraction_coeff(0,2,ωs2)
    diagrams = get_diagrams(DiagramNode((3,1)))
    ω_test= split_freqs_into_bubbles(ωs, diagrams[2])
    corr = diagram_correction(ω_test)
    print(corr)
    corr.poly
    corr.exponent == contraction_coeff(2,0,ωs2).exponents[1]
    corr.exponent  ∈ [1 + ω_d^2, ω_0]
    corr + corr
    ω_1 == ω_1
    isequal(-0.0, ω_1 - ω_1)
    simplify(1.0/ω_1 + -1/ω_1)

    
    k=2
    @syms t::Real
    gs, Vs = simplify(effective_hamiltonian(k, ωs2, hs2))
    @show Vs
    @show gs
    @show [typeof(g*V) for (g,V) in zip(gs, Vs)]
    @show [g*V for (g,V) in zip(gs, Vs)]
    @show (typeof(a'*a)==QuantumCumulants.QMul{Nothing})
    sums = 0
    for (g,V) in zip(gs, Vs)
        if (typeof(g*V) != QuantumCumulants.QMul{Nothing})
            sums+= 0
        else
            sums += g*V
        end
    end
    @show(simplify(sums))
    (typeof(1) == QuantumCumulants.QMul{Nothing}) ? print("yes") : print("no")
    #collect_rule = @rule(+(~~xs) => ~~xs)
    #H_list = simplify(collect_rule(H))
    #typeof(ω_d*σz1)
    
    secondOrderRule = @acrule((~a)*(~y)*(~z) + (~b)*(~y)*(~z) => ((~a)+(~b))*(~y)*(~z))
    firstOrderRule = @acrule((~a)*(~y) + (~b)*(~y) => ((~a)+(~b))*(~y))
    rules = Chain([firstOrderRule,secondOrderRule])
    H1 = Fixpoint(rules)(sums)
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
    

    l = [1,2,3]
    w = [4,5,6]
    # problem with 0 frequencies and limits
"""
end