using QuantumCumulants
using QuantumGraining

ha = NLevelSpace(Symbol(:atom));     hc = FockSpace(:cavity);  
h = hc ⊗ hc                                                     # Atom and Field spaces

σ(i,j) = Transition(ha, Symbol("σ_{$k}"), i, j)    
σz = 2*σ(1,1) - 1
σx = σ(1,2) + σ(2,1)                                            # Spin operators

@cnumbers ω_c ω_a Ω

ωs = [ω_c+ω_a, -ω_c-ω_a, -ω_c+ω_a, ω_c-ω_a]
hs = Ω/2 .* [a*σ(2,1), a'*σ(1,2), a*σ(1,2), a'*σ(2,1)]

H, L_list, J_list = effective_lindblad(h_list, ω_list, order=2)



