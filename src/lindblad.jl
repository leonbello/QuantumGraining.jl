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

"""
    repeated_combinations(h::Vector, g::Vector, Ω::Vector, n::Int)

Generate all possible combinations of elements from the input vectors `h`, `g`, and `Ω` with repetition, forming combinations of length `n`.

# Arguments
- `h::Vector`: Vector of elements for combination from `h`
- `g::Vector`: Vector of elements for combination from `g`
- `Ω::Vector`: Vector of elements for combination from `Ω`
- `n::Int`: Length of combinations to generate

# Returns
- `perm_h::Vector`: Vector of combinations of elements from `h`
- `perm_g::Vector`: Vector of combinations of elements from `g`
- `perm_Ω::Vector`: Vector of combinations of elements from `Ω`
"""
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

"""
    effective_hamiltonian(h::Vector, g::Vector{Number}, Ω::Vector{Number}, k::Int)
    
Given a truncation order `k`, a list of frequencies `Ω`, couplings `g`, and operators `h` representing the raw Hamiltonian, 
this function returns new frequencies, coupling strengths, and operators representing the Hamiltonian to order `k`.

## Arguments
- `h::Vector`: A vector of operators representing the raw Hamiltonian.
- `g::Vector{Number}`: A vector of coupling strengths.
- `Ω::Vector{Number}`: A vector of frequencies.
- `k::Int`: The truncation order.

## Returns
- `ops_eff::Vector`: A vector of operators representing the effective Hamiltonian to order `k`.
- `merged_gs::Vector`: A vector of merged coupling strengths.
- `ωs_eff::Vector`: A vector of frequencies representing the effective Hamiltonian to order `k`.

## Details
The function calculates the effective Hamiltonian by considering all possible combinations of operators, frequencies, and coupling strengths up to order `k`. 
It simplifies the expressions and merges duplicate coupling strengths.

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
end    

"""
    effective_hamiltonian(h::Vector, gs::Vector, Ω::Vector, k::Int; as_dict=false, remove_constants=true)

Compute the effective Hamiltonian for a given system.

## Arguments
- `h::Vector`: A vector of Hamiltonian terms.
- `gs::Vector`: A vector of coupling strengths.
- `Ω::Vector`: A vector of frequencies.
- `k::Int`: The number of terms to consider.

## Keyword Arguments
- `as_dict::Bool=false`: If `true`, the output will be returned as a dictionary of operators and frequencies.
- `remove_constants::Bool=true`: If `true`, remove constant terms from the output.

## Returns
- If `as_dict` is `true`, returns a tuple `(gs_eff, ωs_eff)` where `gs_eff` is a dictionary of operators and `ωs_eff` is a dictionary of frequencies.
- If `as_dict` is `false`, returns a tuple `(unique_hs, unique_gs, ωs_eff)` where `unique_hs` is a vector of unique operators, `unique_gs` is a vector of unique coupling strengths, and `ωs_eff` is a vector of frequencies.

"""
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

"""
    effective_dissipator_term(h::Vector, gs::Vector, Ω::Vector, k::Int)

Compute the effective dissipator term for a given set of parameters.

# Arguments
- `h::Vector`: Vector of coefficients for the Hamiltonian terms.
- `gs::Vector`: Vector of coefficients for the jump operators.
- `Ω::Vector`: Vector of coefficients for the frequencies.
- `k::Int`: Number of terms in the Hamiltonian.

# Returns
- `J_list::Vector`: Vector of tuples representing the jump operators.
- `γ_list::Vector`: Vector of complex numbers representing the dissipation coefficients.
- `ω_list::Vector`: Vector of sums of frequencies.

"""
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

"""
    effective_dissipator(h::Vector, gs::Vector, Ω::Vector, k::Int; as_dict=true)

Compute the effective dissipator for a given set of parameters.

## Arguments
- `h::Vector`: A vector of operators representing the raw Hamiltonian.
- `gs::Vector`: A vector of coupling strengths.
- `Ω::Vector`: A vector of frequencies.
- `k::Int`: The truncation order.
- `as_dict::Bool`: (optional) If `true`, returns the dissipator as a dictionary of operators and their corresponding dissipator terms. If `false`, returns the dissipator as separate vectors for operators, dissipator terms, and frequencies. Default is `true`.

## Returns
- If `as_dict` is `true`, returns a tuple `(γs_dict, ωs_dict)` where `γs_dict` is a dictionary of operators and their corresponding dissipator terms, and `ωs_dict` is a dictionary of operators and their corresponding frequencies.
- If `as_dict` is `false`, returns a tuple `(ops_eff, γs_eff, ωs_eff)` where `ops_eff` is a vector of operators representing the effective dissipator, `γs_eff` is a vector of dissipator terms, and `ωs_eff` is a vector of frequencies.

"""
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



"""
    drop_high_freqs(freqs_list::Vector, freqs_subs, cutoff=0.1)

Given a list of frequencies and a list of substitutions, returns only the low frequencies.

# Arguments
- `freqs_list::Vector`: A list of frequencies.
- `freqs_subs`: A list of substitutions.
- `cutoff=0.1`: The cutoff value for determining low frequencies.

# Returns
- `rwa`: A list of indices corresponding to the low frequencies.
- `freqs_low`: A list of low frequencies.
"""
function drop_high_freqs(freqs_list::Vector, freqs_subs, cutoff=0.1)
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


"""
    drop_high_freqs(gs_dict, freqs_dict, freqs_vals; cutoff=0.1)

Given dictionaries holding the frequencies and couplings, drops all high-frequency contributions.

## Arguments
- `gs_dict`: A dictionary holding the coupling strengths.
- `freqs_dict`: A dictionary holding the frequencies.
- `freqs_vals`: A dictionary of substitutions for the frequencies.
- `cutoff=0.1`: The cutoff value for determining high frequencies.

## Returns
- `gs_dict`: The updated dictionary of coupling strengths after dropping high-frequency contributions.
- `freqs_dict`: The updated dictionary of frequencies after dropping high-frequency contributions.
"""
function drop_high_freqs(gs_dict, freqs_dict, freqs_vals; cutoff=0.1)
    for (key, freq) in freqs_dict
        num_val = substitute(freq, freqs_vals)
        if !isapprox(round(abs(num_val.val)), 0; rtol=cutoff)
            delete!(freqs_dict, key)
            delete!(gs_dict, key)
        end
    end
    return gs_dict, freqs_dict
end

"""
    gaussian_to_cutoff(gs, freq_vals; cutoff=0.1, keep_small_exponents=true)

Apply a cutoff to Gaussian exponents based on their absolute value.

# Arguments
- `gs`: A list of `ContractionCoefficient` objects representing Gaussian functions.
- `freq_vals`: A dictionary mapping frequency keys to their corresponding values.
- `cutoff`: The cutoff value for the absolute value of the Gaussian exponents. Default is `0.1`.
- `keep_small_exponents`: A boolean indicating whether to keep small exponents or set them to zero. Default is `true`.

# Returns
- A list of `ContractionCoefficient` objects with the cutoff applied to the Gaussian exponents.

"""
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

"""
    gaussian_to_cutoff(gs::Dict, ωs::Dict, freq_vals; cutoff=0.1, keep_small_exponents=true)

Apply a cutoff to Gaussian exponents based on their absolute value.

# Arguments
- `gs`: A dictionary mapping frequency keys to `ContractionCoefficient` objects representing Gaussian functions.
- `ωs`: A dictionary mapping frequency keys to their corresponding values.
- `freq_vals`: A dictionary mapping frequency keys to their corresponding values.
- `cutoff`: The cutoff value for the absolute value of the Gaussian exponents. Default is `0.1`.
- `keep_small_exponents`: A boolean indicating whether to keep small exponents or set them to zero. Default is `true`.

# Returns
- A dictionary mapping frequency keys to `ContractionCoefficient` objects with the cutoff applied to the Gaussian exponents.

"""
function gaussian_to_cutoff(gs::Dict, ωs::Dict, freq_vals; cutoff=0.1, keep_small_exponents=true)
    gs_low = Dict(key => gs[key] for key in collect(keys(ωs)) if haskey(gs, key))
    gs_list = gaussian_to_cutoff(collect(values(gs_low)), freq_vals; 
                                    cutoff=cutoff, keep_small_exponents=keep_small_exponents)
    gs_low = Dict(key => gs_list[i] for (i, key) in enumerate(collect(keys(ωs))))

    return gs_low
end