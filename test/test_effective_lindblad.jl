using Revise
using Symbolics
using Test
using QuantumGraining
using QuantumCumulants

#@testset "corrections" begin
module Tst
    using Test
    using QuantumGraining 
    using QuantumCumulants   
    using Symbolics
    
    # Setup
    h_cav = FockSpace(:cavity)
    h_atom = NLevelSpace(:atom, (:g,:e))
    h = tensor(h_cav, h_atom)

    @qnumbers a::Destroy(h) σ::Transition(h)

    σz = 2* σ(:e, :e)  - 1
    σx = σ(:g, :e) + σ(:e, :g)
    σy = 1im* (σ(:e, :g) - σ(:g, :e))    

    @variables ω_c ω_a 

    ωs_rabi = [ω_c + ω_a, -ω_c - ω_a, -ω_c + ω_a, ω_c - ω_a]
    hs_rabi = [a*σ(:e, :g),a'*σ(:g, :e),a*σ(:g, :e), a'*σ(:e, :g)]

    γ_list, ω_list, J_list, Jd_list = effective_dissipator(hs_rabi, ωs_rabi, 2)
    size(γ_list)
    size(ω_list)
    size(J_list)
    size(Jd_list)

    unique(J_list)
    unique(Jd_list)
end