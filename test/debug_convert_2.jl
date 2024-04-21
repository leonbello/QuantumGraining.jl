using Revise
using QuantumCumulants
using QuantumOptics
using QuantumGraining
using Plots
using Symbolics
using LaTeXStrings
using IterTools


### QuantumCumulants.jl definitions
# Define hilbert space
hf_qc = FockSpace(:cavity)
ha_qc = NLevelSpace(:atom,(:g,:e))
h_qc = hf_qc ⊗ ha_qc

# Define the fundamental operators and couplings
@qnumbers a_qc::Destroy(h_qc) σ_qc::Transition(h_qc)
σp_qc = σ_qc(:e, :g)
σm_qc = σ_qc(:g, :e)

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

Id = [I_c, I_a]
hs_qc = [a_qc, a_qc', σp_qc, σm_qc]
hs_qo = [a_qo, a_qo', σp_qo, σm_qo]

fieldnames(typeof(2*a_qc*σm_qc))
fieldnames(typeof(a_qc))

a_qc.aon
σm_qc.aon
fieldnames(typeof(a_qc.hilbert))
a_qc.hilbert.spaces
σp_qc.hilbert.spaces

# function qnumber_to_qop(qn::QuantumCumulants.QMul, op_subs, Id; mul = tensor)
#     num_spaces = length(qn.args_nc[1].hilbert.spaces)
#     fac = qn.arg_c
#     spaces = [[] for i in 1:num_spaces]
#     for op in qn.args_nc
#         push!(spaces[op.aon], op)
#     end

#     new_ops = []
#     for ops in spaces
#         new_ops = push!(new_ops, [substitute(op, op_subs) for op in ops])
#     end
    
#     @show new_ops
#     for i in eachindex(new_ops)
#         if !isempty(new_ops[i])
#             new_ops[i] = fac*prod(new_ops[i])
#         else
#             new_ops[i] = Id[i]
#         end
#     end

#     return mul(new_ops...)
# end

# fieldnames(typeof(a_qc + a_qc))
# function qnumber_to_qop(qn::QuantumCumulants.QAdd, op_subs, Id; mul = tensor)
#     qops = [qnumber_to_qop(qm, op_subs, Id; mul = mul) for qm in qn.arguments]
#     return sum(qops)
# end

# function qnumber_to_qop(qn::QuantumCumulants.QSym, op_subs, Id; mul = tensor)
#     return 1//2*qnumber_to_qop(2*qn, op_subs, Id; mul = mul)
# end


σm_qc isa QuantumCumulants.QSym
2*σm_qc isa QuantumCumulants.QTerm

test1 = qnumber_to_qop(2*a_qc'*a_qc*σm_qc, hs_qc .=> hs_qo, Id)
test2 = qnumber_to_qop(a_qc'*a_qc, hs_qc .=> hs_qo, Id)
test3 = qnumber_to_qop(2*a_qc'*a_qc*σm_qc + a_qc'*a_qc + a_qc, hs_qc .=> hs_qo, Id)
test4 = qnumber_to_qop(a_qc, hs_qc .=> hs_qo, Id)

