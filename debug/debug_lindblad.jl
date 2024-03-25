using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining

@variables g ωc ωa
Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
gvec = (g/2).*[1, 1, 1, 1]

# Hilbert space definitions
h_cav = FockSpace(:cavity)
h_atom = NLevelSpace(:atom, (:g,:e))
h = tensor(h_cav, h_atom)

# Operator definitions
@qnumbers a::Destroy(h) σ::Transition(h)
σm = σ(:e, :g)
σp = σ(:g, :e)
hvec = [a*σm, a'*σp, a*σp, a'*σm]

ops_eff, Ω_eff, g_eff = effective_hamiltonian_term(hvec, gvec, Ω, 1)

ops_eff, Ω_eff, g_eff = effective_hamiltonian(hvec, gvec, Ω, 2)



### RWA
freqs_subs = Dict(
    ωa => 1,
    ωc => 1.01
)
rwa, Ω_rwa = drop_high_freqs(Ω_eff, freqs_subs)
g_rwa = gaussian_to_cutoff(g_eff[rwa], freqs_subs, keep_small_exponents=true)
ops_rwa = ops_eff[rwa]

@variables t τ
H_rwa_2 = sum(symbolic_hamiltonian(g_rwa, ops_rwa, Ω_rwa, t, τ))

@cnumbers g_qc ωc_qc ωa_qc τ_qc
@syms t::Real

subs = Dict(
    g => g_qc,
    ωc => ωc_qc,
    ωa => ωa_qc
)
H_rwa_qc = qc_convert(g_rwa, ops_rwa, Ω_rwa, subs, t, τ_qc)