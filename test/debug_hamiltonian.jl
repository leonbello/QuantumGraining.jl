using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining

@variables g ωc ωa
Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
gvec = (g/2).*[1, 1, 1, 1]

h_cav = FockSpace(:cavity)
h_atom = NLevelSpace(:atom, (:g,:e))
h = tensor(h_cav, h_atom)

@qnumbers a::Destroy(h) σ::Transition(h)
σm = σ(:e, :g)
σp = σ(:g, :e)
σz = σ(:e, :e)
hvec = [a*σm, a'*σp, a*σp, a'*σm]

Σ = ωa + ωc
Δ = ωc - ωa

ops_eff, g_eff, Ω_eff = effective_hamiltonian_term(hvec, gvec, Ω, 2)  

@show a*a'
@show σm*σp
@show σp*σm

@show a'*σm*a*σp
@show a*σp*a'*σm
@show a'*σp*a*σm
@show a*σm*a'*σp

@show ops_eff
unique_hs, unique_gs, unique_ωs = expand_operators(ops_eff, g_eff, Ω_eff)


ids = findall(x -> isequal(x, a'*a), unique_hs)
c_test = simplify_contraction(sum(unique_gs[ids]))
@show c_test.polys

gs, ωs = group_operators(unique_hs, unique_gs, unique_ωs)

@show unique_hs

@test issetequal(gs[a'*a].exponents, c_test.exponents)
@test issetequal(gs[a'*a].prefacs, c_test.prefacs)

g_eff, Ω_eff = effective_hamiltonian(hvec, gvec, Ω, 2; as_dict=true)
@show g_eff

@show g_eff[a'*a].exponents
@show g_eff[a'*a].prefacs

@show g_eff[σ(:e,:e)].exponents
@show g_eff[σ(:e,:e)].prefacs

@show g_eff[a'*a*σ(:e,:e)].exponents
@show g_eff[a'*a*σ(:e,:e)].prefacs

g1_test = g^2/4*contraction_coeff(2, 0, -[-ωa - ωc, ωa + ωc])
g2_test = g^2/4*contraction_coeff(2, 0, -[ωa - ωc, ωc - ωa])
@show g1_test
@show g2_test

@show (g2_test - g1_test).prefacs
@show g_eff[a'*a].prefacs

@show -(g2_test - g1_test).prefacs
@show g_eff[σ(:e,:e)].prefacs

@show 2*(g2_test + g1_test).prefacs
@show g_eff[a'*a*σ(:e,:e)].prefacs

for op in unique_hs
    if op ≠ 1
        println(op)
        @test issetequal(g_eff[op].exponents, gs[op].exponents)
        @test issetequal(simplify.(g_eff[op].prefacs .- gs[op].prefacs), [0,0])
    end
end

# group_operators and simplify_contraction, as well as expand_operators seems to work
# The operator-ordered don't match what I would expet from theory.
# Specifically, the a'*a contributions should cancel, but in the output the still seem to exist.

# function I may have screwed up (not limited to): 
    # +(ContractionCoefficient, ContractionCoefficient)
    # merge_duplicate_exponents()
    # group_unique_operators() - for some reason the ContractionCoefficient is corrupted 