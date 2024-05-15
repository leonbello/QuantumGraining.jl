using Symbolics
using QuantumCumulants
using QuantumOptics

function contraction_to_function(g, ω, ps)
    @variables t, τ
    expr = sum(to_symbol(g, τ).*exp(1im*ω*t))
    func_re = build_function(expr.re, [ps..., τ, t]...; expression=false)
    func_im = build_function(expr.im, [ps..., τ, t]...; expression=false)
    
    func(args, t) = func_re(args..., t) + 1im*func_im(args..., t)
    return func
end

function lindblad_function(gs, Ωs, γs, ωs, h_src, h_tgt, ps)
    n = length(h_src)
    Id = h_tgt[(n+1):end]
    h_tgt = h_tgt[1:n]
    
    
    @variables t, τ
    op_subs = h_src .=> h_tgt
    J = []
    Jdagger = []
    rate_funcs = []
    for (J_val, γ) in γs
        J1, J2 = J_val

        push!(J, qnumber_to_qop(J1, op_subs, Id))
        push!(Jdagger, qnumber_to_qop(J2, op_subs, Id))    
        push!(rate_funcs, contraction_to_function(γ, ωs[J_val], ps))
    end
    
    function L_func(t, ρ; args)
        H = hamiltonian_function(gs, Ωs, h_src, [h_tgt..., Id...], ps)
        rates = map(f -> f(args, t), rate_funcs)
        return (H(t, ρ; args=args), J, Jdagger, rates) 
    end

    return L_func
end


function hamiltonian_function(gs, ωs, h_src, h_tgt, ps)
    n = length(h_src)
    Id = h_tgt[(n+1):end]
    h_tgt = h_tgt[1:n]
    
    func_gs = []
    @variables t, τ

    func_gs = [contraction_to_function(g, ω, ps) for (g, ω) in zip(values(gs), values(ωs))]
    
    new_ops = []
    for op in keys(gs)
        push!(new_ops, qnumber_to_qop(op, h_src .=> h_tgt, Id))
    end

    function H_func(t, ψ; args)
        H = []
        for (f, op) in zip(func_gs, new_ops)
            g = f(args, t)
            push!(H, g*op)
        end
        return sum(H)
    end
    return H_func
end

function qnumber_to_qop(qn::QuantumCumulants.QMul, op_subs, Id; mul = tensor)
    num_spaces = length(qn.args_nc[1].hilbert.spaces)
    fac = qn.arg_c
    spaces = [[] for i in 1:num_spaces]
    for op in qn.args_nc
        push!(spaces[op.aon], op)
    end

    new_ops = []
    for ops in spaces
        new_ops = push!(new_ops, [substitute(op, op_subs) for op in ops])
    end
    
    for i in eachindex(new_ops)
        if !isempty(new_ops[i])
            new_ops[i] = prod(new_ops[i])
        else
            new_ops[i] = Id[i]
        end
    end

    return fac*mul(new_ops...)
end

function qnumber_to_qop(qn::QuantumCumulants.QAdd, op_subs, Id; mul = tensor)
    qops = [qnumber_to_qop(qm, op_subs, Id; mul = mul) for qm in qn.arguments]
    return sum(qops)
end

function qnumber_to_qop(qn::QuantumCumulants.QSym, op_subs, Id; mul = tensor)
    return 1//2*qnumber_to_qop(2*qn, op_subs, Id; mul = mul)
end

### Possibly redundant ###
function qc_convert(gs::Vector, ops::Vector, freqs::Vector, cnumbers_dict, t, τ)
    gs_qc = substitute(to_symbol.(gs, τ), cnumbers_dict)
    freqs_qc = [substitute(ω, cnumbers_dict) for ω in freqs]
    ops_qc = [substitute(op, cnumbers_dict) for op in ops]
    
    return [g*exp(-im*ω*t)*op for (g, ω, op) in zip(gs_qc, freqs_qc, ops_qc)]
end

function qc_convert(gs, ωs, cnumbers_dict, t, τ) 
    ops_list = gs_list = collect(keys(gs))
    gs_list = collect(values(gs))
    ωs_list = collect(values(ωs))
    return qc_convert(gs_list, ops_list, ωs_list, cnumbers_dict, t, τ)
end

function normal_ordered_dictionary(h_src, h_tgt; order=1)
    hd_src, hs_src = h_src
    hd_tgt, hs_tgt = h_tgt

    ops_src = []
    ops_tgt = []
    for (l, r) in [(i, order - i) for i in 0:order]
        curr_src = vec(prod.(collect(product(repeat([hd_src], l)..., repeat([hs_src], r)...))))
        curr_tgt = vec(prod.(collect(product(repeat([hd_tgt], l)..., repeat([hs_tgt], r)...))))
        push!(ops_src, curr_src...)
        push!(ops_tgt, curr_tgt...)
    end
    return Dict(ops_src .=> ops_tgt)
end

function convert_expressions(gs, ωs, h_src, h_tgt, p_subs, τ; order=1, as_dict=true)
    p_subs = Dict(p_subs)
    # check expand_operators and group_operators

    ops_subs = normal_ordered_dictionary(h_src, h_tgt; order=order)    

    #new_ops = [substitute(op, ops_subs) for op in old_ops]

    new_gs = []
    new_ωs = []
    new_ops = []
    old_ops = []
    for op in keys(ops_subs)
        if haskey(gs, op)
            push!(new_gs, substitute(to_symbol(gs[op], τ), p_subs).val)
            push!(new_ωs, substitute(ωs[op], p_subs).val)
            push!(new_ops, ops_subs[op])
            push!(old_ops, op)
        end
    end

    gs_dict = Dict(old_ops .=> new_gs)
    ωs_dict = Dict(old_ops .=> new_ωs)
    if as_dict
        ops_dict = Dict(old_ops .=> new_ops)
        return ops_dict, gs_dict, ωs_dict
    else
        return new_ops, new_gs, new_ωs
    end
end