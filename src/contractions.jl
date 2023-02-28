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
    split_freqs(freqs::Array, ububs::Int, dbubs::Int)

Splits a frequency array into an array of `ububs` up-bubbles and `dbubs` down-bubbles 

Arguments:
    `fres`  - unified array of frequencies, should be of the form [μ1, μ2, ..., μl, ν1, ν2, ..., νr]
    `l` - number of up-modes
    `r` - number of down-modes

Returns:
    freqs_up - array of only the up-modes
    freqs_dn - array of only the down-modes
"""
function split_freqs(freqs::Array, l::Int, r::Int)
    freqs_up = freqs[1:l]
    freqs_dn = freqs[l+1:l+r]
    return freqs_up, freqs_dn
end

"""
    count_modes(diagram)
Counts the total number of up-modes and down-modes in a given diagram
"""
function count_modes(diagram)
    bubs = tuple(map(sum, zip(diagram...)))[1]
    return bubs[1], bubs[2] 
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
    #println(freqs_up)
    #println(freqs_dn)
    freqs_dn = reverse(freqs_dn)
    ind_μ = 0
    ind_ν = 0
    for bubble in diagram
        #println(bubble)
        (μ_len, ν_len) = bubble
        push!(μ, freqs_up[ind_μ+1:ind_μ+μ_len])
        ind_μ += μ_len
        push!(ν, freqs_dn[ind_ν+1:ind_ν+ν_len]) 
        ind_ν += ν_len
        #println(μ)
        #println(ν)
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
function contraction_coeff(order::Tuple{Int, Int}, ω::Array)
    left, right = order
    return contraction_coeff(left, right, ω)
end

function contraction_coeff(left::Int, right::Int, freqs::Array)
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)
    c = 0
    c_list = []
    
    for diagram in diagrams
        # reversing since the diagrams is ordered left-to-right instead of right-to-left
        #println(reverse(diagram))
        ω = split_freqs_into_bubbles(freqs, reverse(diagram))
        #ω_old = reverse(split_freqs_into_bubbles(freqs, diagram))
        #println(ω)
        push!(c_list, (diagram_correction(ω), diagram))
        c += c_list[end][1] 
    end
    return c, c_list
end

"""
    diagram_correction(ω)


Gives an expression for the effective diagram correction.

    - `ω`: A list of vector (μi, νi) with the frequency values for each mode.
    # Arguments

    # Returns
    - A symbolic expression for the diagram contribution.
    # Example
    - For a diagram d = [(3, 2), (2, 1)] we would ω = [(μ1, ν1), (μ2, ν2)] where μi and νi are vectors containing mode values.
 """   
function diagram_correction(ω)
    num_bubbles = length(ω)
    (s_list, stag_list) = find_all_poles(ω)                             # singular indices for current bubble
    println(s_list)
    println(stag_list)
    total_num_poles = count_poles(s_list, stag_list)

    sols = find_integer_solutions(3*num_bubbles, total_num_poles)    
    unl_list = reshape_sols(sols, total_num_poles, num_bubbles)         
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
    - `total_num_poles`: the total number of singular poles in the diagram
    - `bubble_idx`: the index of the bubble for which we want to calculate the correction factor.
    - `s`: a list of the singular poles in the upper modes.
    - `stag`: a list of the singular poles in the lower modes.

# Returns 
    - an array of all bubble factors.
"""
function calculate_bubble_factor(ω, bubble_idx, sols, s, stag)
    μ, ν = ω[bubble_idx]   
    ν = reverse(ν)                   
    @cnumbers τ

    f(x) = exp(-0.5*τ^2*x^2)
    # finite part of the bubble factor (not including the expansion terms)
    start_idx = (bubble_idx == 1) ? 2 : 1
    
    sum_μ = isempty(μ) ? 0 : sum(μ)
    sum_ν = isempty(ν) ? 0 : sum(ν)  # explicity deal with the case where one the vectors is empty (up-bubble or down-bubble)
    prefac = -f(sum_μ + sum_ν)/(vec_factorial(μ[end:-1:start_idx], include_poles = false)*vec_factorial(ν, include_poles=false))

    return prefac*sum( singular_expansion( μ[start_idx:end], ν, sols, s, stag, first_bubble = (bubble_idx == 1) ) )
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
        return prod(s)*prod(stag)
    end
end


function singular_expansion(μ, ν, sols, s, stag; first_bubble = false)
    l = length(μ)
    r = length(ν)
    @cnumbers τ
    
    ju_list = filter(x -> !(x in s), 1:l)
    
    jl_list = filter(x -> !(x in stag), 1:r)  # non-singular indices
    terms = []
        # first inner sum
    for (idx, (n, u, d)) in enumerate(sols)
        analytic_terms = []
        
        start_idx = first_bubble ? 2 : 1
        for k in 0:floor(Int, n/2)
            sum_μ = isempty(μ) ? 0 : sum(μ)  
            sum_ν = isempty(ν) ? 0 : sum(ν)

            freq_sum = isequal(sum_μ + sum_ν, 0) && (n - 2*k) == 0 ? 1 : (sum_μ + sum_ν)
            l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases

            push!(analytic_terms, taylor_coeff(n, k)/Float64(factorial(n))*τ^(2*(n - k))*(freq_sum)^(n - 2*k)*(l_plus_r)^n)
        end
        # second inner sum
        # mu_list and ml_list hold vectors of mu values
        mu_list = find_integer_solutions(length(ju_list), u)
        ml_list = find_integer_solutions(length(jl_list), d)

        # denominator normalization factor - equals to 1 if s or s' is empty.
        # does not depend on m and mtag, so we can pull it out of the sum

        denominator = calculate_normalization(s, stag)
        pole_terms = []
        for (mu_vec, ml_vec) in product(mu_list, ml_list)   # for each solution vector   

            fac_u = []
            for ju in ju_list
                idx_u = indexin(ju, ju_list)[1]             # pick up the corresponding index for m_{Ju}

                push!(fac_u, norm_fac(μ, ju - start_idx + 1, mu_vec[idx_u]))
            end
            prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

            fac_v = []
            for jl in jl_list
                idx_l = indexin(jl, jl_list)[1]             # pick up the corresponding index for m_{Jl}
                push!(fac_v, norm_fac(ν, jl, ml_vec[idx_l]))
            end
            prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

            push!(pole_terms, prod_fac_u*prod_fac_v/denominator) 
        end
        poles_sum = isempty(pole_terms) ? 0 : sum(pole_terms)
        analytic_sum = isempty(analytic_terms) ? 0 : sum(analytic_terms)
        push!(terms, analytic_sum*poles_sum)
    end
    return terms
end

#= old code, may be rdundant =#
##-----------------OLD FUNCTIONS--------------------##

#NOTE: contraction_coefficients is used in effective_hamiltonian - not redundant


"""
    _Ω(freqs...)
Simple helper function to calculate the normalization factor for the correction terms.
"""
function _Ω(freqs)
    n = length(freqs)
    terms = [sum(freqs[1:(n - i)]) for i in 0:(n-1)]
    return n > 0 ? prod(terms) : 1
end

"""
    get_max_modes(diagram)

Returns the maximal number of modes in a diagram.
        
# Examples
    diagram = [(3,4), (2, 5)]
    n, m = get_max_modes(diagram)  
    > n = 3
    > m = 5
"""
function get_max_modes(diagram::Array{Tuple{Int64, Int64}})
    ububs = zeros(Int64, length(diagram))
    dbubs = zeros(Int64, length(diagram))
    for (i, bubble) in enumerate(diagram)
        ububs[i] = bubble[1]
        dbubs[i] = bubble[2]
    end
    return (maximum(ububs), maximum(dbubs))
end

"""
    bubble_coeff(l, r, endmode)

    Calculates the coefficient of a single bubble, 
    returns a symbolic function that evaluates to the coefficient as a function of frequency.

    Arguments:
    - μ: a vector of the symbolic up-mode frequencies in the bubble.
    - ν: a vector of the symbolic down-mode frequencies in the bubble.

    Returns:
    - f_{l,r}(μ, ν): a symbolic expression for a single bubble coeffeicient.
"""
function bubble_coeff(l, r, endmode=0) 
    @cnumbers τ
    f(ω, τ) = exp(-1/2*ω^2*τ^2)                      # define filter function

    # define an expression for the bubble correction with the frequencies (μ, ν) as function arguments
    function _c(μ, ν) 
        sum_μ = l > 0 ? sum(μ) : 0
        sum_ν = r > 0 ? sum(ν) : 0
        return -f(sum_μ + sum_ν, τ)/_Ω(μ[1:(end-endmode)])*_Ω(ν)
    end
    return _c                                                                    # should return a symbolic *function*
end
