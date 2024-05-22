using Symbolics

"""
Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
"""

"""
    repeated_combinations(arr::Array, n::Int)
    
Given a collection of terms and a number `n`, returns all possible `n` repetitions of the elements in the collection.
"""
function repeated_combinations(h::Vector, Ω::Vector, n::Int)
    perm_h = []
    perm_Ω = []
    total_combinations = length(h)^n

    # Generate combinations
    for i in 1:total_combinations
        combination_index = digits(i-1; base=length(h), pad=n)    # Convert the integer i to a combination index in base `length(arr)` with `n` digits
        curr_h = []
        curr_Ω = []

        # Use the combination index to pick elements from the input array and form the combination
        for index in combination_index
            push!(curr_h, h[index + 1])  
            push!(curr_Ω, Ω[index + 1]) 
        end

        # Add the current combination to the list of combinations
        push!(perm_Ω, curr_Ω)
        push!(perm_h, curr_h)
    end

    return perm_h, perm_Ω
end

function repeated_combinations(h::Vector, g::Vector, Ω::Vector, n::Int)
    perm_h = []
    perm_Ω = []
    perm_g = []
    total_combinations = length(h)^n

    # Generate combinations
    for i in 1:total_combinations
        combination_index = digits(i-1; base=length(h), pad=n)    # Convert the integer i to a combination index in base `length(arr)` with `n` digits
        curr_h = []
        curr_Ω = []
        curr_g = []

        # Use the combination index to pick elements from the input array and form the combination
        for index in combination_index
            push!(curr_h, h[index + 1])  
            push!(curr_Ω, Ω[index + 1]) 
            push!(curr_g, g[index + 1])  
        end

        # Add the current combination to the list of combinations
        push!(perm_Ω, curr_Ω)
        push!(perm_h, curr_h)
        push!(perm_g, curr_g)
    end

    return perm_h, perm_g, perm_Ω
end

# function symbolic_eff_ham(gs, ops, ωs,T)
#     return to_symbol(gs)*ops*exp(-1im*sum(ωs)*T)
# end

"""
    effective_hamiltonian(h::Vector, g::Vector{Number}, Ω::Vector{Number}, k::Int)
Given a truncation order `k`, and a list of frequencies `Ω`, couplings `g` and operators `h`` representing the raw Hamiltonian, 
returns new frequencies, coupling strengths and operators representing the Hamiltonian to order `k`.
"""
function effective_hamiltonian_term(h::Vector, gs::Vector, Ω::Vector, k::Int)
    perm_h, perm_g, perm_Ω = repeated_combinations(h, gs, Ω, k)
    
    ops_eff = []
    ωs_eff  = []
    gs_eff  = []

    idx = []
    for i in eachindex(perm_h)
        ω = perm_Ω[i]
        if length(ω) > 1 && isequal(ω[1] + ω[2], 0)
            push!(idx, i)
        end
        op = perm_h[i]
        g = simplify_contraction(1//2*prod(perm_g[i])*(contraction_coeff(k, 0, ω) + contraction_coeff(k, 0, -reverse(ω))))
        push!(ops_eff, op)                    
        push!(ωs_eff, ω)
        push!(gs_eff, g)
    end
    if length(ops_eff[1]) == 1
        ops_eff = vcat(ops_eff...)
        ωs_eff = vcat(ωs_eff...)
        gs_eff = vcat(gs_eff...)
    end

    if k > 1
        ops_eff = [simplify(prod(op)) for op in ops_eff]
        ωs_eff = [sum(ω) for ω in ωs_eff]
    end

    merged_gs = [merge_duplicate_exponents(g) for g in gs_eff]
    return ops_eff, merged_gs, ωs_eff
    #return ops_eff, merge_duplicate_exponents.(gs_eff), ωs_eff
end    

function effective_hamiltonian(h::Vector, gs::Vector, Ω::Vector, k::Int; as_dict=false, remove_constants=true)
    ops_eff = []
    ωs_eff  = []
    gs_eff  = []

    for i in 1:k
        op, g, ω = effective_hamiltonian_term(h, gs, Ω, i)
        push!(ops_eff, op...)
        push!(ωs_eff, ω...)
        push!(gs_eff, g...)
    end
    
    unique_hs, unique_gs, unique_ωs = expand_operators(ops_eff, gs_eff, ωs_eff)
    
    if as_dict
        gs_eff, ωs_eff = group_operators(unique_hs, unique_gs, unique_ωs; as_dict=true)
        if remove_constants
            for key in keys(gs_eff)
                if key isa Number 
                    delete!(gs_eff, key)
                    delete!(ωs_eff, key)
                end
            end
        end
        return gs_eff, ωs_eff
    else
        unique_hs, unique_gs = group_operators(unique_hs, unique_gs, unique_ωs; as_dict=false) 
        return unique_hs, unique_gs, ωs_eff
    end

end

"""
* effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
# To be deprecated
function effective_dissipator_term(h::Vector, Ω::Vector, k::Int)
    γ_list  = []
    ω_list  = []
    J_list  = []
    
    for i in 1:k  
        perm_h, perm_Ω = repeated_combinations(h, Ω, i)

        for (ω, h) in zip(perm_Ω, perm_h) 
            for l in 1:(i-1)
                γ = -im*(contraction_coeff(l, k-l, ω) - contraction_coeff(k-l, l, -reverse(ω)))
                
                J = prod(h[1:l]) # Jump operators
                L = prod(h[l+1:k])

                push!(γ_list, γ) 
                push!(ω_list, sum(ω))
                push!(J_list, (J, L))
            end
        end
    end
    return γ_list, ω_list, J_list
end

function effective_dissipator_term(h::Vector, gs::Vector, Ω::Vector, k::Int)
    γ_list  = []
    ω_list  = []
    J_list  = []
    
    for i in 1:k  
        perm_h, perm_g, perm_Ω = repeated_combinations(h, gs, Ω, k)

        for (ω, g, h) in zip(perm_Ω, perm_g, perm_h) 
            for l in 1:(i-1)
                γ = -im*(contraction_coeff(l, k-l, ω) - contraction_coeff(k-l, l, -reverse(ω)))
                γ = merge_duplicate_exponents(γ)
                J = prod(h[1:l]) # Jump operators
                L = prod(h[l+1:k])

                push!(γ_list, prod(g)*γ)
                push!(ω_list, sum(ω))
                push!(J_list, (J, L))
            end
        end
    end
    return J_list, γ_list, ω_list
end

function effective_dissipator(h::Vector, gs::Vector, Ω::Vector, k::Int; as_dict=true)
    ops_eff = []
    ωs_eff  = []
    γs_eff  = []

    for i in 1:k
        op, γ, ω = effective_dissipator_term(h, gs, Ω, i)
        push!(ops_eff, op...)
        push!(ωs_eff, ω...)
        push!(γs_eff, γ...)
    end
    
    #unique_hs, unique_γs, unique_ωs = expand_operators(ops_eff, γs_eff, ωs_eff)

    if as_dict
        #γs_eff, ωs_eff = group_operators(unique_hs, unique_γs, unique_ωs; as_dict=true)
        keys = []
        values_γ = []
        values_ω = []
        for i in eachindex(ops_eff)
            if !issetequal(γs_eff[i].prefacs, zeros(length(γs_eff[i].prefacs)))
                push!(keys, ops_eff[i])
                push!(values_γ, γs_eff[i])
                push!(values_ω, ωs_eff[i])
            end
        end
        γs_dict = Dict(keys .=> values_γ)
        ωs_dict = Dict(keys .=> values_ω)
        #γs_dict = Dict(ops_eff[i] => γs_eff[i] for i in eachindex(ops_eff))
        #ωs_dict = Dict(ops_eff[i] => ωs_eff[i] for i in eachindex(ops_eff))
        return γs_dict, ωs_dict
    else
        return ops_Eff, γs_eff, ωs_eff
    end
end

function expand_operators(hs)
    """
        Goes over any sum of operators and breaks it into the constituent operators while preserving the order of the other vectors
    """
    unique_hs = []

    for h in hs
        new_h, _ = expand_operator(h)
        push!(unique_hs, new_h...)
    end

    return unique_hs
end

function expand_operators(hs, gs, ωs)
    """
        Goes over any sum of operators and breaks it into the constituent operators while preserving the order of the other vectors
    """
    unique_hs = []
    unique_gs = []
    unique_ωs = []

    for (h, g, ω) in zip(hs, gs, ωs)
        new_h, facs = expand_operator(h)

        push!(unique_hs, new_h...)
        push!(unique_gs, [fac*g for fac in facs]...)
        push!(unique_ωs, ω*ones(length(new_h))...)
    end

    return unique_hs, unique_gs, unique_ωs
end

function expand_operator(h)
    ops = []
    facs = []
    if h ≠ 0
        if h isa QuantumCumulants.QAdd  # if it's a sum of operators
            for op in h.arguments
                push!(ops, (op isa QuantumCumulants.QMul) ? prod(op.args_nc) : op)
                push!(facs, (op isa QuantumCumulants.QMul) ? op.arg_c : 1)
            end
        elseif h isa QuantumCumulants.QMul
            push!(ops, h)
            push!(facs, 1)
        end
    end
    return ops, facs
end

function group_operators(hs, gs, ωs; as_dict=true)
    new_hs = []
    new_gs = []
    new_ωs = []

    for i in eachindex(hs)
        idx = findfirst(x -> isequal(x, hs[i]), new_hs)
        if !isnothing(idx)
            # If it's a duplicate, find its index in new_hs i=14
            new_gs[idx] += gs[i]
        else
            # If it's not a duplicate, add hs[i] to new_hs
            push!(new_hs, hs[i])
            # Add gs[i] to new_gs and add corresponding element from ws to new_ws
            push!(new_gs, gs[i])
            push!(new_ωs, ωs[i])
        end
    end

    new_gs = simplify_contraction.(new_gs)

    if as_dict
        gs_dict = Dict(new_hs[i] => new_gs[i] for i in eachindex(new_hs))
        ωs_dict = Dict(new_hs[i] => new_ωs[i] for i in eachindex(new_hs))
        return gs_dict, ωs_dict
    else
        return new_hs, new_gs, new_ωs
    end
end

function group_operators(hs)
    new_hs = []

    for i in eachindex(hs)
        idx = findfirst(x -> isequal(x, hs[i]), new_hs)
        if isnothing(idx)
            push!(new_hs, hs[i])
        end
    end

    return new_hs

end

"""
* effective_lindblad(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_lindblad(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_lindblad(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
# function effective_lindblad(k::Int, ω::Array, h::Array)
#     return effective_hamiltonian(k, ω, h), effective_dissipator(k, ω, h)
# end

function drop_high_freqs(freqs_list::Vector, freqs_subs, cutoff=0.1)
    """
        Given a list of frequencies, and a list of substitutions, returns only the low frequencies.
    """
    freqs_low = []
    rwa = []
    for (i, freq) in enumerate(freqs_list)
        freq_val = substitute(freq, freqs_subs)
        if isapprox(round(abs(freq_val.val)), 0; rtol=cutoff)    
            push!(rwa, i)
            push!(freqs_low, freq)
        end
    end
    return rwa, freqs_low
end

function drop_high_freqs(gs_dict, freqs_dict, freqs_vals; cutoff=0.1)
    """
        Given dictionaries holding the frequencies and couplings,
        drops all high-frequency contributions.
    """

    for (key, freq) in freqs_dict
        num_val = substitute(freq, freqs_vals)
        if !isapprox(round(abs(num_val.val)), 0; rtol=cutoff)
            delete!(freqs_dict, key)
            delete!(gs_dict, key)
        end
    end
    return gs_dict, freqs_dict
end

function symbolic_hamiltonian(gs::Vector, ops::Vector, Ω::Vector, t, τ)
    terms = []
    for (g, op, ω) in zip(gs, ops, Ω)
        if isequal(ω, 0)
            ft = 1
        else
            ft = Symbolics.Term(exp, [im*ω*t])
        end

        term = to_symbol(g, τ)*ft
        if !isequal(term, 0)
            push!(terms, (term)*op)
        end
    end
    return terms
end

function symbolic_hamiltonian(gs::Dict, Ω::Dict, t, τ)
    gs_list = collect(values(gs))
    ops_list = collect(keys(gs))
    Ω_list = collect(values(Ω))

    return symbolic_hamiltonian(gs_list, ops_list, Ω_list, t, τ)
end

function gaussian_to_cutoff(gs, freq_vals; cutoff=0.1, keep_small_exponents=true)
    new_gs = []
    for g in gs
        new_exponents = []
        include_ids = []
        for (i, exponent) in enumerate(g.exponents)
            exponent_val = substitute(exponent, freq_vals)
            if isapprox(round(abs(exponent_val.val)), 0; rtol=cutoff)
                push!(new_exponents, keep_small_exponents ? exponent : 0)
                push!(include_ids, i)
            end
        end
        push!(new_gs, ContractionCoefficient(new_exponents, g.prefacs[include_ids], g.polys[include_ids]))
    end
    return merge_duplicate_exponents.(new_gs)
end

function gaussian_to_cutoff(gs::Dict, ωs::Dict, freq_vals; cutoff=0.1, keep_small_exponents=true)
    gs_low = Dict(key => gs[key] for key in collect(keys(ωs)) if haskey(gs, key))
    gs_list = gaussian_to_cutoff(collect(values(gs_low)), freq_vals; 
                                    cutoff=cutoff, keep_small_exponents=keep_small_exponents)
    gs_low = Dict(key => gs_list[i] for (i, key) in enumerate(collect(keys(ωs))))

    return gs_low
end