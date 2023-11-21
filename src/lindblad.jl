using Symbolics
#"""
#Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
#"""
"""
Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
"""

"""
    repeated_combinations(arr::Array, n::Int)
    
Given a collection of terms and a number `n`, returns all possible `n` repetitions of the elements in the collection.
"""
# Helper function (NOTE: not sure what the protocol is to steal something from StackOverflow)
# repeated_combinations(arr::Vector, n::Int) = 
#     [ getindex.(Ref(arr), 1 .+ digits(i-1; base=length(arr), pad=n)) for i=1:length(arr)^n ]
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
            # Add the element at the selected index to the current combination
            push!(curr_h, h[index + 1])  
            push!(curr_Ω, Ω[index + 1]) # Add 1 to the index because Julia arrays are 1-based
        end

        # Add the current combination to the list of combinations
        push!(perm_Ω, curr_Ω)
        push!(perm_h, curr_h)
    end

    return perm_h, perm_Ω
end

"""
    normal_order(ops::Vector)
Given an array of operators, returns the normal ordered product of them.
"""
function normal_order(ops::Vector)
    # if we use QuantumCumulants, the product is enough since they would be automatically normal ordered.
    return prod(ops)
end



function symbolic_eff_ham(gs, ops, ωs,T)
    return to_symbol(gs)*ops*exp(-1im*sum(ωs)*T)
end

"""
    effective_hamiltonian(h::Vector, g::Vector{Number}, Ω::Vector{Number}, k::Int)
Given a truncation order `k`, and a list of frequencies `Ω`, couplings `g` and operators `h`` representing the raw Hamiltonian, 
returns new frequencies, coupling strengths and operators represneting the new Hamiltonian.
"""
# NOTE: This outputs term of the kth order only
function effective_hamiltonian(h::Vector, Ω::Vector, k::Int)
    perm_h, perm_Ω = repeated_combinations(h, Ω, k)
    
    ops_eff = []
    ωs_eff  = []
    gs_eff  = []

    for (i, perm) in enumerate(perm_h)      
        push!(ops_eff, perm) # not sure this works, may be better to keep ops_eff have only unique operators
        ω = perm_Ω[i]
        push!(ωs_eff, ω)
        push!(gs_eff, contraction_coeff(k, 0, ω))
    end
    if length(ops_eff[1]) == 1
        ops_eff = vcat(ops_eff...)
        ωs_eff = vcat(ωs_eff...)
    end
    return ops_eff,ωs_eff, gs_eff
end    

"""
* effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
function effective_dissipator(h::Vector, Ω::Vector, k::Int)
    γ_list  = []
    ω_list  = []
    J_list  = []
    Jd_list = []
    
    for i in 1:k  
        perm_h, perm_Ω = repeated_combinations(h, Ω, i)

        for (ω, h) in zip(perm_Ω, perm_h) 
            for l in 1:(i-1)
                γ = (contraction_coeff(l, k-l, ω) - contraction_coeff(k-l, l, ω))
                
                J = prod(h[1:l]) # Jump operators
                Jd = prod(h[l+1:k])

                push!(γ_list, γ) 
                push!(ω_list, sum(ω))
                push!(J_list, J)
                push!(Jd_list, Jd)
            end
        end
    end
    return γ_list, ω_list, J_list, Jd_list
end

# function effective_dissipator(k::Int, ω::Array, h::Array, fmt=:QuantumCumulants)



"""
* effective_lindblad(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_lindblad(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_lindblad(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
function effective_lindblad(k::Int, ω::Array, h::Array)
    return effective_hamiltonian(k, ω, h), effective_dissipator(k, ω, h)
end