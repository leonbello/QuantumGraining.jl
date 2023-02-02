"""
Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
"""


# helper function for the innermost sum
norm_fac(v, j) = return j == 0 ? 0 : (-j/sum(v[1:j]))


"""
    diagram_correction(ω)

Gives an expression for the effective diagram correction.

    # Arguments
    - `ω`: A list of vector (μi, νi) with the frequency values for each mode.

    # Returns
    - A symbolic expression for the diagram contribution.

    # Example
    - For a diagram d = [(3, 2), (2, 1)] we would ω = [(μ1, ν1), (μ2, ν2)] where μi and νi are vectors containing mode values.
    
"""
function diagram_correction(ω)
    # find all singularities
    num_bubbles = length(ω)
    (s_list, stag_list) = find_all_poles(ω)                             # singular indices for current bubble
    total_num_poles = count_poles(s_list, stag_list)

    sols = find_integer_solutions(3*num_bubbles, total_num_poles)
    unl_list = sols     
    #unl_list = reshape_sols(sols, total_num_poles, num_bubbles)         # partition for the inner sum -- ##problem here##

    bubble_factors = []                                                 # array holding the terms of the outer sum
    l_tot = 0
    for (i, (μ, ν)) in enumerate(ω)
        fac = calculate_bubble_factor(ω, i, unl_list[:, i], s_list[i], stag_list[i])
        push!(bubble_factors, fac)
        l_tot += length(μ)
    end
    return (-1)^l_tot*prod(bubble_factors)
end



"""
    calculate_bubble_factor(ω, bubble_idx, total_num_poles, s, stag)
Returns the bubble factor for a single bubble.

# Arguments
    - `ω`: the frequency values for the whole diagram.
    - `bubble_idx`: the index of the bubble for which we want to calculate the correction factor.
    - `total_num_poles`: the total number of singular poles in the diagram
    - `s`: a list of the singular poles in the upper modes.
    - `stag`: a list of the singular poles in the lower modes.

# Returns 
    - an array of all bubble factors.
"""
function calculate_bubble_factor(ω, bubble_idx, sols, s, stag)
    μ, ν = ω[bubble_idx]                      
    @cnumbers τ
    f(x) = exp(-1/2*τ^2*x^2)

    # finite part of the bubble factor (not including the expansion terms)
    
    start_idx = (bubble_idx == 1) ? 2 : 1
    sum_μ = isempty(μ) ? 0 : sum(μ)
    sum_ν = isempty(ν) ? 0 : sum(ν)  # explicity deal with the case where one the vectors is empty (up-bubble or down-bubble)
    prefac = -f(sum_μ + sum_ν)/(vec_factorial(μ[end:-1:start_idx], include_poles = false)*vec_factorial(ν, include_poles=false))

    return prefac*sum(singular_expansion(μ, ν, sols, s, stag))
end

# Helper function for singular expansion
function calculate_normalization(s, stag)
    if isempty(s) && isempty(stag)
        return 1
    elseif isempty(s) && !isempty(stag)
        return prod(stag)
    elseif isempty(stag) && !isempty(s)
        return prod(s)
    elseif !isempty(s) && !isempty(stag)
        return prod([su*sl for (su, sl) in (s, stag)])
    end
end

function singular_expansion(μ, ν, sols, s, stag)
    l = length(μ)
    r = length(ν)
    @cnumbers τ
    
    ju_list = filter(x -> !(x in s), 1:l)
    jl_list = filter(x -> !(x in stag), 1:r)  # non-singular indices
    
    terms = []
    for (idx, (n, u, d)) in enumerate(sols)
        # first inner sum
        analytic_terms = []
        for k in 0:floor(Int, n/2)
            sum_μ = isempty(μ) ? 0 : sum(μ)
            sum_ν = isempty(ν) ? 0 : sum(ν)
            freq_sum = isequal(sum_μ + sum_ν, 0) && (n - 2*k) == 0 ? 1 : (sum_μ + sum_ν)
            l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases

            push!(analytic_terms, taylor_coeff(n, k)/factorial(n)*τ^(2*(n - k))*(freq_sum)^(n - 2*k)*(l_plus_r)^n)
        end
        # second inner sum
        # mu_list and ml_list hold vectors of mu values
        mu_list = find_integer_solutions(length(ju_list), u)
        ml_list = find_integer_solutions(length(jl_list), d)

        # denominator normalization factor - equals to 1 if s or s' is empty.
        denominator = calculate_normalization(s, stag)

        pole_terms = []
        for (mu_vec, ml_vec) in product(mu_list, ml_list)   # for each solution vector   
            numerator = []
            for (ju, jl) in (ju_list, jl_list)              # for each non-singular factor
                idx_u = indexin(ju, ju_list)[1]             # pick up the corresponding index for m_{Ju}
                idx_l = indexin(jl, jl_list)[1]
                fac = (norm_fac(μ, ju))^mu_vec[idx_u]*(norm_fac(ν, jl))^ml_vec[idx_l]
                # make sure the first mode in the first bubble is omitted
                push!(numerator, fac)
            end
            prod_numerator = isempty(numerator) ? 1 : prod(numerator)
            push!(pole_terms, prod_numerator/denominator) 
        end
        poles_sum = sum(pole_terms)
        push!(terms, sum(analytic_terms)*poles_sum)
    end
    return terms
end


"""
    repeated_combinations(arr::Array, n::Int)
    
        Given a collection of terms and a number `n`, returns all possible `n` repetitions of the elements in the collection.
"""
# Helper function (NOTE: not sure what the protocol is to steal something from StackOverflow)
repeated_combinations(arr::Vector, n::Int) = 
    [ getindex.(Ref(arr), 1 .+ digits(i-1; base=length(arr), pad=n)) for i=1:length(arr)^n ]


"""
    effective_hamiltonian(k::Int, ω::Array, h::Array)
    
        Given a Hamiltonian as a list of operators and corresponding frequencies,
        returns the effective TCG Hamiltonian up to order `k`.
"""
function effective_hamiltonian(k::Int, ω::Array, h::Array)
    # Functionality to sum over all contractions of up to kth order 
    # w and h are arrays of symbolic cnumbers from QuantumCumulants
    g_list = []
    V_list = []
    h_eff = 0
    ω_list = repeated_combinations(ω, k)
    h_list = repeated_combinations(h, k)
    for i in 1:k
        g = 0
        V = 0
        lc(ω) = (contraction_coeff(k, 0))(ω)                # returns a function of μ and ν
        rc(ω) = (contraction_coeff(0, k))(ω)                
        for (ω, h) in zip(ω_list, h_list) 
            g += 1//2*(lc(ω) + rc(ω))
            V += prod(h) #? 
        end
        push!(g_list, g) 
        push!(V_list, V)
        h_eff += g*V
    end
    return h_eff
    # effective_ham = 0
    # for n in 1:k
    #     ω_list = repeated_combinations(ω,n)
    #     h_list = repeated_combinations(h,n)
    #     effective_ham += sum([effective_hamiltonian(diagram,ω_list,h_list) for diagram in get_diagrams(DiagramNode((n,0)))])    
    # end
    # return effective_ham
end


# """
# * effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
# * effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
# * effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
# """

# function effective_dissipator(k::Int, ω::Array, h::Array, fmt=:QuantumCumulants)



"""
* effective_lindblad(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_lindblad(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_lindblad(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""
