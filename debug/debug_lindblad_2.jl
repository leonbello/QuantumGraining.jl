using Symbolics
using Revise
using QuantumCumulants
using QuantumOptics
using QuantumGraining
using OrdinaryDiffEq
using ModelingToolkit
using Plots

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


# First-order TCG (RWA)
g_eff_1, Ω_eff_1 = effective_hamiltonian(hs, gs, Ω, 1; as_dict=true)

freqs_subs = Dict(
    ωa => 1,
    ωc => 1.01
)
Ω_low_1 = drop_high_freqs(Ω_eff_1, freqs_subs)
g_low_1 = gaussian_to_cutoff(g_eff_1, Ω_low_1, freqs_subs; keep_small_exponents=false)

@cnumbers g_qc ωc_qc ωa_qc τ_qc
@syms t_qc::Real

subs = Dict(
    g => g_qc,
    ωc => ωc_qc,
    ωa => ωa_qc
)
H_qc_1 = qc_convert(g_low_1, Ω_low_1, subs, t_qc, τ_qc)
@show sum(H_qc_1)

# Second-order TCG
hs_eff, gs_eff, ωs_eff = effective_hamiltonian_term(hs, gs, Ω, 2)

g_eff_2, Ω_eff_2 = effective_hamiltonian(hs, gs, Ω, 2; as_dict=true)

@show g_eff_2[a'*a]

freqs_subs = Dict(
    ωa => 1,
    ωc => 1.01
)
Ω_low_2 = drop_high_freqs(Ω_eff_2, freqs_subs)
g_low_2 = gaussian_to_cutoff(g_eff_2, Ω_low_2, freqs_subs; keep_small_exponents=true)

expand_operator(a - 3*a' + 2*a'*a*σ(:e,:e))
expand_operators([a + a', a - a'], [1, 1], [3, 5])

J_eff_2, ωs_eff_2 = effective_dissipator(hs, gs, Ω, 2)

@show J_eff_2

H_tcg_qc_1 = qc_convert(g_low_1, Ω_low_1, subs, t, τ_qc)


# convert_expression
function convert_expression(gs, ωs, ops_dict, params_dict)
    new_ops = [substitute(op, ops_dict) for op in keys(gs) ]
    new_gs = [substitute(to_symbol.(g, τ), params_dict) for g in values(gs)]
    new_ωs = [substitute(to_symbol.(ω, τ), params_dict) for ω in values(ωs)]

    gs_dict = Dict(new_ops .=> new_gs)
    ωs_dict = Dict(new_ops .=>new_ωs)

    return gs_dict, ωs_dict
end

