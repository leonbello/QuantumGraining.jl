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

    for i in eachindex(perm_h)
        ω = perm_Ω[i]      
        push!(ops_eff, perm_h[i])                    
        push!(ωs_eff, ω)
        push!(gs_eff, prod(1//2*perm_g[i])*(contraction_coeff(k, 0, ω) + contraction_coeff(k, 0, -reverse(ω))))
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

    return ops_eff, merge_duplicate_exponents.(gs_eff), ωs_eff
end    

function effective_hamiltonian(h::Vector, gs::Vector, Ω::Vector, k::Int)
    ops_eff = []
    ωs_eff  = []
    gs_eff  = []

    for i in 1:k
        op, g, ω = effective_hamiltonian_term(h, gs, Ω, i)
        push!(ops_eff, op...)
        push!(ωs_eff, ω...)
        push!(gs_eff, g...)
    end

    return ops_eff, gs_eff, ωs_eff
end

"""
* effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
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



"""
* effective_lindblad(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_lindblad(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_lindblad(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
# function effective_lindblad(k::Int, ω::Array, h::Array)
#     return effective_hamiltonian(k, ω, h), effective_dissipator(k, ω, h)
# end

function drop_high_freqs(freqs_list, freqs_subs, cutoff=0.1)
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

function symbolic_hamiltonian(gs, ops, Ω, t, τ)
    terms = []
    for (g, op, ω) in zip(gs, ops, Ω)
        if isequal(ω, 0)
            ft = 1
        else
            ft = Symbolics.Term(exp, [-im*ω*t])
        end

        term = to_symbol(g, τ)*ft
        println(term)
        if !isequal(term, 0)
            push!(terms, term*op)
        end
    end
    return terms
end

function gaussian_to_cutoff(gs, freq_vals; cutoff=0.1, keep_small_exponents=false)
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

function qc_convert(gs, ops, freqs, cnumbers_dict, t, τ)
    gs_qc = substitute(to_symbol.(gs, τ), cnumbers_dict)
    freqs_qc = [substitute(ω, cnumbers_dict) for ω in freqs]
    ops_qc = [substitute(op, cnumbers_dict) for op in ops]
    
    return sum([g*exp(-im*ω*t)*op for (g, ω, op) in zip(gs_qc, freqs_qc, ops_qc)])
end