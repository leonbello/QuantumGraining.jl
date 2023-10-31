#Test file for lindblad.jl 
using QuantumCumulants
using IterTools
using Symbolics
using SymbolicUtils.Rewriters
using SymPy
using SymbolicUtils
using QuantumGraining

module Tst
    using Test
    
    #Setup-
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
    
    @syms t::Real
    eff_ham_rab, ops_eff_rab, ωs_eff_rab, gs_eff_rab = effective_hamiltonian(hs_rab, ωs_rab, 2, t)
    sum(eff_ham_rab)
    ops_eff_rab
    render(latexify(to_symbol(gs_eff_rab[14])))
    
    #Numerical test
    ω_n = 1
    ω_m = 1.0001

    
    ωs_rab_num = [ω_n + ω_m, -ω_n - ω_m, -ω_n + ω_m, ω_n - ω_m]
    hs_rab_num = [a*σ(:e, :g),a'*σ(:g, :e),a*σ(:g, :e), a'*σ(:e, :g)]
    
    @syms t::Real
    eff_ham_rab_num, ops_eff_rab_num, ωs_eff_rab_num, gs_eff_rab_num = effective_hamiltonian(hs_rab_num, ωs_rab_num, 2, t)

    gs_eff_rab_num[2].polys
    render(latexify(to_symbol(gs_eff_rab_num[2])))
    eff_ham_rab_num

    to_symbol(contraction_coeff(2,0,[ω_a,ω_c]))
    to_symbol(contraction_coeff(0,2,[ω_a,ω_c]))
end