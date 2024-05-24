
## Short example
As a simple illustrative example, let us consider the implementation of a second-order time coarse graining of the Rabi model:
```
using QuantumCumulants
using QuantumGraining

@variables g ωc ωa
Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
gvec = (g/2).*[1, 1, 1, 1]

# Hilbert space definitions (QuantumCumulants.jl)
h_cav = FockSpace(:cavity)
h_atom = NLevelSpace(:atom, (:g,:e))
h = tensor(h_cav, h_atom)

# Operator definitions
@qnumbers a::Destroy(h) σ::Transition(h)
σm = σ(:e, :g)
σp = σ(:g, :e)
hvec = [a*σm, a'*σp, a*σp, a'*σm]

order=2
g_eff, Ω_eff = drop_high_freqs(effective_hamiltonian(hvec, gvec, Ω, order; as_dict=true)..., Dict(ωa => 1, ωc => 1.01))
γ_eff, ω_eff = drop_high_freqs(effective_dissipator(hvec, gvec, Ω, order)..., Dict(ωa => 1, ωc => 1.01)) 
```

The code above returns an effective Lindbladian that generates the time-coarse grained evolution of the Rabi-model up to second-order. The input is a list of frequencies and their corresponding operators, and the output is a new list of operators, frequencies and coupling strengths. The `drop_high_freqs` function is used to remove the high-frequency terms from the effective Hamiltonian and dissipator, to simplify the resulting expressions.
