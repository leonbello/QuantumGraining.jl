# QuantumGraining

**QuantumGraining.jl** offers a practical, generalized approach to the time-coarse graining (TCG) method in quantum optics. Often in quantum optics, we are not interested in the short time-scale dynamics -- they are either trivial, or unmeasurable, and keeping track of them is very computationally expensive. The time-coarse graning approach gives an effective description for the time-coarse grained dynamics, while preserving the slow time-scale effects of the short time-scale dynamics.

The main challenge with the TCG approach is that the calculations are very involved, with the number of terms growing rapidly with the order of truncation. **QuantumGraining.jl** automates this process, representing the effective Lindbladians in terms of abstract operators that are easily integrated into other symbolic packages, such as **QuantumCumulants.jl** and **QuantumOptics.jl**.

* The model (Hamiltonian) is specified, along with the time-coarse graining time-scale.
* The effective Lindbladian is calculated, using an efficient, recursive approach. The resulting Lindbladian is stored stored as a symbolic expression using the [**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) framework, which is also used for any additional simplification and rewriting.
* Finally, the symbolic Hamiltonian can be solved in **QuantumOptics.jl** or using **QuantumCumulants.jl**. 

# Development Status
Note that **QuantumGraining.jl** is still at an early stage of development.
[![Build Status](https://github.com/leonbello/QuantumGraining.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/leonbello/QuantumGraining.jl/actions/workflows/CI.yml?query=branch%3Amain)

