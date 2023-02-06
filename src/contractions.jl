using QuantumCumulants
using SymbolicUtils


# helper function for the innermost sum
norm_fac(v, j) = return j == 0 ? 0 : (-j/sum(v[1:j]))


"""
    contraction_coeff(left::Int, right::Int)

    Calculates the coefficient of a whole contraction, given the contraction and input frequencies

    Arguments:
    - left: the left-order of the contraction
    - right: the right-order of the contraction
    - freqs: the frequency array

    Returns:
    - c: a symbolic expression for the contraction coeffeicient.
"""

function contraction_coeff(left::Int, right::Int, freqs::Array)
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)           #Get all the diagrams for the contraction (left, right)
    c = 0
    for diagram in diagrams
        μ, ν = [],[]
        ind = 0
        for (i, bubble) in enumerate(diagram)
            (μ_len, ν_len) = bubble
            push!(μ, freqs[ind+1:ind+μ_len])
            ind += μ_len
            push!(ν,freqs[ind+1:ind+ν_len]) 
            ind += ν_len
        end
        ω = tuple.(μ, ν)

        c += diagram_correction(ω)
    end
    return c
end




function contraction_coeff(order::Tuple{Int, Int}, ω::Array)
    left, right = order
    return contraction_coeff(left, right, ω)
end


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
    unl_list = reshape_sols(sols, total_num_poles, num_bubbles)           # partition for the inner sum

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
        return prod([su*sl for (su, sl) in product(s, stag)])
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
            for (ju, jl) in product(ju_list, jl_list)       # for each non-singular factor
                idx_u = indexin(ju, ju_list)[1]             # pick up the corresponding index for m_{Ju}
                idx_l = indexin(jl, jl_list)[1]
                fac = (norm_fac(μ, ju))^mu_vec[idx_u]*(norm_fac(ν, jl))^ml_vec[idx_l]
                push!(numerator, fac)
            end
            prod_numerator = isempty(numerator) ? 1 : prod(numerator)
            push!(pole_terms, prod_numerator/denominator) # don't need an array for the pole terms
        end
        poles_sum = isempty(pole_terms) ? 0 : sum(pole_terms)
        push!(terms, sum(analytic_terms)*poles_sum)
    end
    return terms
end


##-----------------OLD FUNCTIONS--------------------##
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


#function contraction_coeff(left::Int, right::Int, divide_freqs = false)
#    node = DiagramNode((left, right))
#    diagrams = get_diagrams(node)
#    function _c(μ, ν)
#        ex = 0                                                                  # initialize an empty expression
#        for diagram in diagrams
#           #Check if the diagram has any poles
#           s_list, stag_list = find_all_poles(ω)
#           total_num_poles = sum([length(su) + length(sl) for (su, sl) in (s_list, stag_list)])
#             for (i, bubble) in enumerate(diagram)
#                    endmode = (i == length(diagram)) ? 1 : 0
#                    diag_ex = diag_ex*(bubble_coeff(left, right, endmode))(μ, ν)    # calculate the bubble correction term and multiply it in
#                end
#
#            else #We calculate the finite term of the Taylor expansion
#                @cnumbers τ
#                f(x) = exp(-1/2*τ^2*x^2)
#                ω = diagram
#                # first sum (n+u+l)_i = total_num_poles
#                unl_list = find_integer_solutions(3*length(ω), total_num_poles)
#                for (i, (μ, ν)) in enumerate(ω)
#                    l = length(μ)
#                    r = length(ν)
#                    
#                    su_list = find_poles(μ)
#                    sv_list = find_poles(ν)                         # singular indices

#                    ju_list = filter(x -> !(x in su_list), 1:l)
#                    jl_list = filter(x -> !(x in sv_list), 1:r)     # non-singular indices
#
#                    num_poles = length(su_list) + length(sv_list)   # total number of poles
#                
#                    # finite part of the bubble factor (not including the expansion terms)
#                    A = -f(sum(μ) + sum(ν))/(vec_factorial(μ, include_poles = false)*vec_factorial(ν, include_poles=false))
#                
#                    unl_list = reshape_sols_vec(find_integer_solutions(3, num_poles), 3, num_poles, length(ω))
#                    
#                    for (idx, (n, u, d)) in enumerate(unl_list[i, :])
#                        numerator_terms = [taylor_coeff(n, k)/factorial(n)*τ^(2*(n - k))*(sum(μ) + sum(ν))^(n - 2*k)*(d + r)^n for k = 1:floor(Int, n/2)]

#                        mu_list = find_integer_solutions(length(ju_list), u)
#                        ml_list = find_integer_solutions(length(jl_list), d)

                        ###
#                        for (mu, ml) in product(mu_list, ml_list)
#                            for (ju, jl) in product(ju_list, jl_list)
#                                non_poles_product = (-ju/sum(μ(1:ju)))^mu_list(ju)*(-jl/sum(ν(1:jl)))^ml_list(jl)
#                            end
                            
#                            poles_product = prod[su*sl for (su, sl) in product(su_list, sl_list)]
#                            denominator_terms = non_poles_product/poles_product
#                        end
                        ###
#                        diagram_factors[idx] = A*sum(numerator_terms)*sum(denominator_terms)
#                    end
#                end
#                diag_ex = prod(diagram_factors)
#            end
#                ex = ex + diag_ex                                                   # add the diagram to the contraction coefficient expression
#        end
#        return ex
#    end
#    if divide_freqs
#        _c(ω) = (-1)^(left)*_c(ω[1:left], ω[left+1:end])  
#    end
#    return _c
#end




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

# """
#     Given a diagram, returns a vector of the terms with the coeffcients for the contraction expression.
# """
# function calculate_coeff(μ, ν, τ, diagram::Array{Tuple{Int64, Int64}})
#     h = []
#     endmode = 0
#     for (i, bubble) in enumerate(diagram)
#         if i == length(diagram)
#             endmode = 1  
#         end
#         n, m = bubble
#         push!(h, coeff(μ[1:n], ν[1:m], τ, endmode))
#     end
#     return prod(h), h
# end

# """
#     Helper function for calculating the coefficient of a single bubble, returns a symbolic expression that evaluates to the coefficient
# """
# function coeff(n, m, endmode=0) 
#     f(ω, τ) = (2*π*τ^2)^(-1/2)*exp(-1/2*ω^2*τ^2)                             # define filter function
#     return f(sum(μ) + sum(ν), τ)/(_Ω(μ)*_Ω(ν)) 
# end
