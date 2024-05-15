using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining
using QuantumOptics

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

order = 3
#ops_eff, Ω_eff, g_eff = effective_hamiltonian_term(hvec, gvec, Ω, order)

g_eff, Ω_eff = effective_hamiltonian(hvec, gvec, Ω, order, as_dict=true)
g_eff, Ω_eff = drop_high_freqs(g_eff, Ω_eff, Dict(ωa => 1, ωc => 1.01))

g_eff

### QuantumOptics.jl definitions
ha_qo = SpinBasis(1//2)
hc_qo = FockBasis(100)
h_qo = hc_qo ⊗ ha_qo

# Operator definitions
σp_qo = sigmap(ha_qo)
σm_qo = sigmam(ha_qo)
a_qo = destroy(hc_qo)
I_a = identityoperator(ha_qo)
I_c = identityoperator(hc_qo)

p_sym = [g, ωc, ωa]

base_qc = [a, a', σm, σp, σ(:e, :e)]
Id = [I_c, I_a]
base_qo = [a_qo, a_qo', σm_qo, σp_qo, σp_qo*σm_qo, Id...]

H_func = hamiltonian_function(g_eff, Ω_eff, base_qc, base_qo, p_sym)

## QuantumOptics.jl
# Units
μs = 1
MHz = 1/μs

tspan = [0:0.01:120μs;]
ψ0 = coherentstate(hc_qo, 4.5) ⊗ spinup(ha_qo)



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