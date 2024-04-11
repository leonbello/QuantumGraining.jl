using Symbolics
using QuantumCumulants
using QuantumGraining

@variables t τ 

@variables g ωc ωa
Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
gs = (g/2).*[1, 1, 1, 1]

# Hilbert space definitions
h_cav = FockSpace(:cavity)
h_atom = NLevelSpace(:atom, (:g,:e))
h = tensor(h_cav, h_atom)

# Operator definitions
@qnumbers a::Destroy(h) σ::Transition(h)
σm = σ(:e, :g)
σp = σ(:g, :e)
hs = [a*σm, a'*σp, a*σp, a'*σm]

hs_eff, gs_eff, ωs_eff = effective_hamiltonian_term(hs, gs, Ω, 2)

g_eff_2, Ω_eff_2 = effective_hamiltonian(hs, gs, Ω, 2; as_dict=true)

@show g_eff_2[a'*a]

freqs_subs = Dict(
    ωa => 1,
    ωc => 1.01
)
Ω_low_2 = drop_high_freqs(Ω_eff_2, freqs_subs)
g_low_2 = gaussian_to_cutoff(g_eff_2, Ω_low_2, freqs_subs; keep_small_exponents=true)

@show simplify_contraction(g_eff_2[a'*a])


####

contraction_coeff(2, 0, [ωc, ωc])
contraction_coeff(2, 0, [-ωc, -ωc])

g = simplify_contraction(contraction_coeff(2, 0, [ωa, ωa]) + contraction_coeff(2, 0, -[ωa, ωa]))
