#Test file for lindblad.jl 
using QuantumCumulants


using Test

diagram = [(2,1), (1,0)]
include("../src/diagrams.jl")
include("../src/contractions.jl")
include("../src/lindblad.jl")
include("../src/printing.jl")

N = 2 #Number of atoms 
h = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k)

σz1 = 2*σ(1,1,1) - 1
σz2 = 2*σ(1,1,2) - 1

σx1 = σ(1,2,1) + σ(2,1,1)
σx2 = σ(1,2,2) + σ(2,1,2)

@cnumbers ω1 ω2 ωd g ϵd

#Define ordered lists of frequencies and hamiltonian terms 
ω_list = [0,0,0,ωd,-ωd]
h_list = [ω1*σz1/2,ω2*σz2/2, g* σx1 * σx2, ϵd* σx2,ϵd* σx2]
@show effective_hamiltonian(2, ω_list, h_list)
