# QuantumGraining

![](docs/imgs/quantumgraining.png)

## Introduction
**QuantumGraining.jl** offers a practical, generalized approach to the time-coarse graining (TCG) method in quantum optics. Often in quantum optics, we are not interested in the short time-scale dynamics -- these are either trivial, or unmeasurable, and keeping track of them is very computationally expensive. The time-coarse graning approach gives an effective, low-frequency description of the time-coarse grained dynamics, without neglecting the effects the high-frequency modes have on this dynamics.

One can think of the time-coarse graining method as a generalization of the rotating-wave approximation, where the time-averaging is performed at the propagator level rather than the Hamiltonian level. More technically, while the rotating-wave approximation neglects all contributions due to counter-propagating terms, the time-coarse graining method is able to include `k`-order terms that can be comprised off these counter-rotating terms,
up to some arbitrary order `k`.  A familiar example would be virtual processes such as the simultaneous emission of a photon together with the excitation of an atom, followed the immediate absorption and de-excitation of the system.

Moreover, the TCG method provides this description in terms of a generalized Linbdladian -- a Hamiltonian and a list of "dissipators" -- that fully generate these processes up to the desired order of approximation. In addition to the obvious benefits of a low-frequency description, this allows one to understand the effects of high-frequency modes directly as emergent dynamics on the relevant modes of interest.

## QuantumGraining.jl

[![Build Status](https://github.com/leonbello/QuantumGraining.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/leonbello/QuantumGraining.jl/actions/workflows/CI.yml?query=branch%3Amain)

The STCG method is very computationally involved, with the number of terms growing rapidly with the order of truncation. The main object of interest is the "contraction coefficient" which encapsulates the relative weight of each possible `k`-order process. Each of these coefficients requires the recursive calculation of lower-order contraction coefficients to compute, a process which is the repeated for each frequency mode in the original Hamiltonian.

<img src="docs/imgs/decomp-tree.png" width=50% height=50%>

**QuantumGraining.jl** automates the TCG process, symbolically calculating the effective Lindbladians in terms of abstract operators that are easily integrated into other symbolic packages, such as **QuantumCumulants.jl** and **QuantumOptics.jl**.

* The model (Hamiltonian) is specified, along with the time-coarse graining time-scale.
* The effective Lindbladian is calculated, using an efficient, recursive approach. The resulting Lindbladian is stored stored as a symbolic expression using the [**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) framework, which is also used for any additional simplification and rewriting.
* Finally, the symbolic Hamiltonian can be solved in **QuantumOptics.jl** or using **QuantumCumulants.jl**. 

## Installation
`QuantumGraining.jl` is  registered in the Julia package registrator, and can be easily installed on using the latest version of Julia. 
1. Within the Julia REPL, press the `]` key. This will switch you to the Pkg REPL, indicated by a prompt like `(@v1.8) pkg>`.
2. In the Pkg REPL, 
    ```
        add QuantumGraining
    ```
3. Since the operator algebra is handled by `QuantumCumulants`, you will need to install it as well,
   ```
        add QuantumCumulants
   ```
5. Once installed, include `QuantumGraining` in your code using:
    ```
        using QuantumGraining
        using QuantumCumulants
    ```
After following these steps, you should have successfully installed `QuantumGraining.jl` and its dependencies.

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
