using QuantumCumulants
using SymbolicUtils

"""
    norm_fac(v, j, mj)
helper function for the innermost sum
"""
function norm_fac(v, j, mj) 
    if j == 0 || mj == 0
        return 1
    else 
        return (-j/sum(v[1:j]))^mj
    end
end 

"""
    split_freqs_into_bubbles(freqs, diagram)
Splits an array of frequencies into an array of tuples of frequencies, matching the dimensions of each bubble in `diagram`.

Argument:
    freqs - an array of frequencies in the form ω = [μ1..., μ2..., ..., μl..., ν1..., ν2..., ..., νr]
    diagram - dimensions of each bubble, where each tuple is the dimension for the corresponding bubble.
Returns:
    ω - an array of frequencies in the form [(μ1, ν1), (μ2, ν2), ..., (μl, νl), ([], ν(l+1)), ..., ([], νr)]
"""
function split_freqs_into_bubbles(freqs, diagram)
    μ, ν = [],[]
    ububs, dbubs = count_modes(diagram)
    freqs_up, freqs_dn = split_freqs(freqs, ububs, dbubs)
    freqs_dn = reverse(freqs_dn)
    ind_μ = 0
    ind_ν = 0
    for bubble in diagram
        (μ_len, ν_len) = bubble
        push!(μ, freqs_up[ind_μ+1:ind_μ+μ_len])
        ind_μ += μ_len
        push!(ν, freqs_dn[ind_ν+1:ind_ν+ν_len]) 
        ind_ν += ν_len
    end
    ω = tuple.(μ, ν)
    return ω
end


"""
    contraction_coeff(left::Int, right::Int)

    Calculates the coefficient of a whole contraction, given the contraction and input frequencies

    Arguments:
    - left: the left-order of the contraction
    - right: the right-order of the contraction
    - freqs: array of frequencies to put in each mode

    Returns:
    - c: a symbolic expression for the contraction coeffeicient.
"""
function contraction_coeff(left::Int, right::Int, freqs::Array)
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)
    c = 0
    c_list = []
    
    for diagram in diagrams
        # reversing since the diagrams is ordered left-to-right instead of right-to-left
        # println(reverse(diagram))
        ω = split_freqs_into_bubbles(freqs, reverse(diagram))
        # ω_old = reverse(split_freqs_into_bubbles(freqs, diagram))
        # println(ω)
        push!(c_list, (diagram_correction(ω), diagram))
        c += c_list[end][1] 
    end
    return c, c_list
end
contraction_coeff(order::Tuple{Int, Int}, ω::Array) = contraction_coeff(order[1], order[2], ω)

"""
    diagram_correction(ω)
Gives an expression for the effective diagram correction.
    
    # Arguments
    - `ω`: A list of vectors (μi, νi) with the frequency values for each mode.

    # Returns
    - A symbolic expression for the diagram contribution.
    # Example
    - For a diagram d = [(3, 2), (2, 1)] we would ω = [(μ1, ν1), (μ2, ν2)] where μi and νi are vectors containing mode values.
 """   
function diagram_correction(d::Diagram{T1, T2}) where {T1, T2}
    num_bubbles = length(d)
    sols = find_integer_solutions(3*num_bubbles, d.num_poles)    
    sols = reshape_sols(sols, d.num_poles, num_bubbles)         
    bubble_factors = []                                                 # array holding the terms of the outer sum
    
    l_tot, r_tot = 0, 0
    for (i, b) in enumerate(d)
        fac = calculate_bubble_factor(b, sols[:, i], d.up_poles[i], d.down_poles[i])
        push!(bubble_factors, fac)
        l_tot += length(b.up)
        r_tot += length(b.down)
    end
    return (-1)^l_tot*prod(bubble_factors)
end
function diagram_correction(ω::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
    d = Diagram(ω)
    return diagram_correction(d)
end


"""
    calculate_bubble_factor(ω, bubble_idx, total_num_poles, s, stag)
Returns the bubble factor for a single bubble.

# Arguments
    - `ω`: the frequency values for the whole diagram.
    - `total_num_poles`: the total number of singular poles in the diagram
    - `bubble_idx`: the index of the bubble for which we want to calculate the correction factor.
    - `s`: a list of the singular poles in the upper modes of the bubble.
    - `stag`: a list of the singular poles in the lower modes of the bubble.

# Returns 
    - an array of all bubble factors.
"""
function calculate_bubble_factor(b::Bubble{T1, T2}, sols, s, stag) where {T1, T2}
    μ, ν = b.up, b.down
    l, r = b.shape
    
    @cnumbers τ
    f(x) = exp(-0.5*τ^2*x^2) # Gaussian filter function

    prefac = -f(sum(μ) + sum(ν))/(vec_factorial(μ)*vec_factorial(ν))
    singular_terms = singular_expansion(μ, ν, sols, s, stag)
    return prefac*sum( singular_terms )
end

function calculate_bubble_factor(ω::Vector{Tuple{Vector{T1}, Vector{T2}}}, bubble_idx::Int, sols, s, stag) where {T1, T2}
    μ, ν = ω[bubble_idx]   
    ν = reverse(ν)
    l = length(μ)
    r = length(ν)
    
    @cnumbers τ
    f(x) = exp(-0.5*τ^2*x^2)
    
    # finite part of the bubble factor (not including the expansion terms)
    first_bubble = (bubble_idx == 1)
    μ0, ν0 = set_indices(l, r, first_bubble)

    sum_μ = isempty(μ) ? 0 : sum(μ)
    sum_ν = isempty(ν) ? 0 : sum(ν)  # explicity deal with the case where one the vectors is empty (up-bubble or down-bubble)
    prefac = -f(sum_μ + sum_ν)/(vec_factorial(μ[end:-1:μ0], include_poles = false)*vec_factorial(ν[end:-1:ν0], include_poles=false))

    return prefac*sum( singular_expansion( μ[μ0:end], ν[ν0:end], sols, s, stag, first_bubble = first_bubble) )
end

function calc_analytic_terms(μ, ν, n)
    l = length(μ)
    r = length(ν)
    @cnumbers τ
    
    analytic_terms = []
    for k in 0:floor(Int, n/2)
        sum_μ = isempty(μ) ? 0 : sum(μ)  
        sum_ν = isempty(ν) ? 0 : sum(ν)
        
        freq_sum = isequal(sum_μ + sum_ν, 0) && (n - 2*k) == 0 ? 1 : (sum_μ + sum_ν)
        l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases

        push!(analytic_terms, taylor_coeff(n, k)/Float64(factorial(n))*τ^(2*(n - k))*(freq_sum)^(n - 2*k)*(l_plus_r)^n)
    end
    return analytic_terms
end

# probably can be deprecated
function set_indices(l, r, first_bubble)
    if first_bubble && !iszero(l)
        μ0, ν0 = 2, 1
    elseif first_bubble && !iszero(r)
        μ0, ν0 = 1, 2
    else
        μ0, ν0 = 1, 1
    end 
end

function calc_pole_terms(μ, ν, s, stag, mu_list, ml_list, ju_list, jl_list, first_bubble)    
    pole_terms = []
    l = length(μ)
    r = length(ν)

    μ0, ν0 = set_indices(l, r, first_bubble)

    # denominator normalization factor - equals to 1 if s or s' is empty.
    # does not depend on m and mtag, so we can pull it out of the sum
    #denominator = calculate_normalization(s, stag)
    denominator = begin
        if isempty(s) && isempty(stag)
            return 1
        elseif isempty(s) && !isempty(stag)
            return prod(stag)
        elseif isempty(stag) && !isempty(s)
            return prod(s)
        elseif !isempty(s) && !isempty(stag)
            return prod(s)*prod(stag)
        end
    end

    for (mu_vec, ml_vec) in product(mu_list, ml_list)   # for each solution vector   

        fac_u = []
        for ju in ju_list
            idx_u = indexin(ju, ju_list)[1]             # pick up the corresponding index for m_{Ju}

            push!(fac_u, norm_fac(μ, ju - μ0 + 1, mu_vec[idx_u]))
        end
        prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

        fac_v = []
        for jl in jl_list
            idx_l = indexin(jl, jl_list)[1]             # pick up the corresponding index for m_{Jl}
            push!(fac_v, norm_fac(ν, jl - ν0 + 1, ml_vec[idx_l]))
        end
        prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

        push!(pole_terms, prod_fac_u*prod_fac_v/denominator) 
    end
    return pole_terms
end

function singular_expansion(μ, ν, sols, s, stag; first_bubble = false)
    l = length(μ)
    r = length(ν)
    @cnumbers τ

    ju_list = filter(x -> !(x in s), 1:l)
    jl_list = filter(x -> !(x in stag), 1:r)  # non-singular indices

    terms = []
    for (idx, (n, u, d)) in enumerate(sols)
        analytic_terms = calc_analytic_terms(μ, ν, n)                                 # first inner sum -- includes all analytic contributions
        
        # mu_list and ml_list hold vectors of mu values
        mu_list = find_integer_solutions(length(ju_list), u)
        ml_list = find_integer_solutions(length(jl_list), d)
        
        pole_terms = calc_pole_terms(μ, ν, s, stag, mu_list, ml_list, ju_list, jl_list, first_bubble)     # second inner sum -- terms due to finite pole contributions     
        
        poles_sum = isempty(pole_terms) ? 0 : sum(pole_terms)
        analytic_sum = isempty(analytic_terms) ? 0 : sum(analytic_terms)
        push!(terms, analytic_sum*poles_sum)
    end
    return terms
end