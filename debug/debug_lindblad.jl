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

order = 4
g_eff, Ω_eff = effective_hamiltonian(hvec, gvec, Ω, order, as_dict=true);
g_low, Ω_low = drop_high_freqs(g_eff, Ω_eff, Dict(ωa => 1, ωc => 1.01));
#g_eff, Ω_eff = drop_high_freqs(g_eff, Ω_eff, Dict(ωa => 1, ωc => 1.01))
ops_list, g_list, Ω_list = effective_hamiltonian(hvec, gvec, Ω, order, as_dict=false);

println(ops_list[9])
@show g_list[9] 
@show simplify(g_list[9].polys[2][3], expand=true, simplify_fractions=false)
val = substitute(g_list[9].polys[2][3], [g, ωc, ωa] .=> [0.2, 2, 2.1])

typeof(g_list[9].polys[2][3].val)
fieldnames(typeof(g_list[9].polys[2][3].val))
@show val

# Units
μs = 1
MHz = 1/μs
args = [2π*0.2MHz, 2π*2MHz, 2π*2.1MHz, 0.2μs]  # g, ωc, ωa, τ 


println("----wo/ operator ordering:----")
println("Prefactors:")
for prefac in g_list[9].prefacs
    val = substitute(prefac, [g, ωc, ωa] .=> [0.2, 2, 2.1])
    @show val
end

println("Exponentials:")
for expon in g_list[9].exponents
    val = substitute(expon, [g, ωc, ωa] .=> [0.2, 2, 2.1])
    @show val
end

println("Polynomials:")
for poly in g_list[9].polys
    println("Polynomial terms:")
    for term in poly
        val = substitute(term, [g, ωc, ωa] .=> [0.2, 2, 2.1])
        @show val
    end
end

println("----w/ operator ordering:----")
println("Prefactors:")
for prefac in g_eff[σ(:e,:e)].prefacs
    val = substitute(prefac, [g, ωc, ωa] .=> [0.2, 2, 2.1])
    @show val
end

println("Exponentials:")
for expon in g_eff[σ(:e,:e)].exponents
    val = substitute(expon, [g, ωc, ωa] .=> [0.2, 2, 2.1])
    @show val
end

println("Polynomials:")
for poly in g_eff[σ(:e,:e)].polys
    println("Polynomial terms:")
    for term in poly
        val = substitute(term, [g, ωc, ωa] .=> [0.2, 2, 2.1])
        @show val
    end
end
