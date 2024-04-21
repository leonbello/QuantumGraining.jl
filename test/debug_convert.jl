using Revise
using QuantumCumulants
using QuantumOptics
using QuantumGraining
using Plots
using Symbolics
using LaTeXStrings
using IterTools

# Parameters
MHz = 1
μs = 1

### QuantumCumulants.jl definitions
# Define hilbert space
hf_qc = FockSpace(:cavity)
ha_qc = NLevelSpace(:atom,(:g,:e))
h_qc = hf_qc ⊗ ha_qc

# Define the fundamental operators and couplings
@qnumbers a_qc::Destroy(h_qc) σ_qc::Transition(h_qc)
σp_qc = σ_qc(:e, :g)
σm_qc = σ_qc(:g, :e)
hs_qc = [a_qc*σm_qc, a_qc'*σp_qc, a_qc*σp_qc, a_qc'*σm_qc]

# Parameters
@syms t_qc::Real
@variables g_qc ωc_qc ωa_qc
p_qc = [g_qc, ωc_qc, ωa_qc]
Ω_qc = [-ωc_qc - ωa_qc, ωc_qc + ωa_qc, -ωc_qc + ωa_qc, ωc_qc - ωa_qc]
gs_qc= (g_qc/2).*[1, 1, 1, 1]

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

hs_qo = [a_qo⊗σm_qo, a_qo'⊗σp_qo, a_qo⊗σp_qo, a_qo'⊗σm_qo]


### TCG
g_eff_1, Ω_eff_1 = effective_hamiltonian(hs_qc, gs_qc, Ω_qc, 1; as_dict=true)
g_eff_2, Ω_eff_2 = effective_hamiltonian(hs_qc, gs_qc, Ω_qc, 2; as_dict=true)

tspan = [0:0.01:80μs;]
ψ0 = coherentstate(hc_qo, 4.5) ⊗ spinup(ha_qo)

base_qc = [a_qc, a_qc', σm_qc, σp_qc, σ_qc(:e, :e)]

Id = [I_c, I_a]
base_qo = [a_qo, a_qo', σm_qo, σp_qo, σp_qo*σm_qo, Id...]


H_1 = hamiltonian_function(g_eff_1, Ω_eff_1, base_qc, base_qo, p_qc)
H_2 = hamiltonian_function(g_eff_2, Ω_eff_2, base_qc, base_qo, p_qc)

args = [2π*0.1MHz, 2π*2MHz, 2π*2.1MHz, 0.2μs]
H_1(0, ψ0; args=args)
H_2(0, ψ0; args=args)
tout1, ψt1 = timeevolution.schroedinger_dynamic(tspan, ψ0, (t, ψ) -> H_1(t, ψ; args=args));
tout2, ψt2 = timeevolution.schroedinger_dynamic(tspan, ψ0, (t, ψ) -> H_2(t, ψ; args=args));

plot(tout1, real(expect(1, a_qo'*a_qo, ψt1)), lw=2.5, label="a'*a - 1");
plot!(tout2, real(expect(1, a_qo'*a_qo, ψt2)), lw=2.5, label="a'*a - 2");

xlabel!("Time [μs]")
ylabel!(L"$\langle n \rangle$")





##


@variables τ
new_ops, new_gs, new_ωs = convert_expressions(g_eff_2, Ω_eff_2, h_src, h_tgt, p_qc .=> p_qc, τ; order=2)
@show new_ops

for (key, val) in new_ops
    println("$key => $(typeof(val)) \n")
end


no_dict = normal_ordered_dictionary((hd_qc, h_qc), (hd_qo, h_qo); order=2)
@show length(values(no_dict))
for e in collect(keys(no_dict))
    if e != 0 && e != 1
        println(e)
    end
end

tcg_operators = [keys(g_eff_1)..., keys(g_eff_2)...]
for e in tcg_operators
    if e != 0 && e != 1
        println(e)
    end
end
issubset(tcg_operators, keys(no_dict))
setdiff(tcg_operators, keys(no_dict))

