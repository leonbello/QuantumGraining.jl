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

"""
    effective_hamiltonian(h::Vector, g::Vector{Number}, Ω::Vector{Number}, k::Int)
Given a truncation order `k`, and a list of frequencies `Ω`, couplings `g` and operators `h`` representing the raw Hamiltonian, 
returns new frequencies, coupling strengths and operators represneting the new Hamiltonian.
"""
# STILL WORK IN PROGRESS, DO NOT RUN
function effective_hamiltonian(h::Vector, Ω::Vector{Number}, k::Int)
    perm_h, perm_Ω = coresponding_combinations(h, Ω, k)
    # [a, a, a] [ω, ω, ω]
    # [a', a, a'] [-ω, ω, -ω]
    h_eff = []
    ωs_eff  = []
    gs_eff  = []

    for i in enumerate(perm_h)        
        push!(h_eff, normal_order(perm_h[i])) # not sure this works, may be better to keep ops_eff have only unique operators
        ω = sum(perm_Ω[i])
        push!(ωs_eff, ω)
        push!(gs_eff, 1/2*(contraction_coeff(k, 0, ω) + contraction_coeff(0, k, ω)))
    end

    return h_eff, ωs_eff, gs_eff
end

"""
    effective_hamiltonian(k::Int, ω::Array, h::Array)
    
        Given a Hamiltonian as a list of operators and corresponding frequencies,
        returns the effective TCG Hamiltonian up to order `k`.
"""
function effective_hamiltonian(K::Int, ω::Array, h::Array)
    # Functionality to sum over all contractions of up to kth order 
    # w and h are arrays of symbolic cnumbers from QuantumCumulants
    g_list = []
    V_list = []
    freq_list = []
    @syms t::Real
    h_eff = 0
    ks = []
    ωs = []
    contractions = []
    for k in 1:K
        ω_list = repeated_combinations(ω, k)
        h_list = repeated_combinations(h, k)
        g = 0
        V = 0               
        for (ω, h) in zip(ω_list, h_list) 
            push!(ks, k)
            push!(ωs, ω)
            g = 0.5*(contraction_coeff((k,0),ω)[1] + contraction_coeff((k,0),-reverse(ω))[1])*exp(1.0im*sum(ω)*t)
            g = simplify(g)
            #println(typeof(g))
            push!(contractions, contraction_coeff((k,0),ω)[1])
            V = prod(h) #? 
            #println(typeof(V))
            push!(g_list, g)
            push!(V_list, V)
            push!(freq_list, sum(ω))
            h_eff += simplify(g*V)
            #(typeof(g*V) == QuantumCumulants.QMul{Nothing}) ? h_eff += g*V : h_eff += 0
        end
    end
    #H_eff = 0
    #terms = tuple.(g_list, V_list)
    #for (j,(g,V)) in enumerate(terms)
    #    @syms f(j,t)
    #    H_eff += g*V*f(j,t)
    #end
    return h_eff, contractions, ks, ωs, g_list, V_list
    # effective_ham = 0
    # for n in 1:k
    #     ω_list = repeated_combinations(ω,n)
    #     h_list = repeated_combinations(h,n)
    #     effective_ham += sum([effective_hamiltonian(diagram,ω_list,h_list) for diagram in get_diagrams(DiagramNode((n,0)))])    
    # end
    # return effective_ham
end


"""
* effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
function effective_dissipator(K::Int, ω::Array, h::Array)
    # Functionality to sum over all contractions of up to kth order 
    # w and h are arrays of symbolic cnumbers from QuantumCumulants
    g_list = []
    J_list = []
    Jd_list = []
    @syms t::Real
    for k in 1:K  
        ω_list = repeated_combinations(ω, k)
        h_list = repeated_combinations(h, k)              
        for (ω, h) in zip(ω_list, h_list) 
            for k1 in 1:(k-1)
                g = (contraction_coeff((k1,k-k1),ω) - contraction_coeff((k-k1,k1),ω))*exp(1im*sum(ω)*t)
                J = prod(h[1:k1]) # Jump operator
                Jd = prod(h[k1+1:k])
                push!(g_list, g) 
                push!(J_list, J)
                push!(Jd_list, Jd)
            end
        end
    end
    return g_list, J_list, Jd_list
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