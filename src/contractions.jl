using QuantumCumulants
using SymbolicUtils
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

"""
    contraction_coeff(left::Int, right::Int)

    Calculates the coefficient of a whole contraction, 
    returns a symbolic function that evaluates to the coefficient as a function of frequency.

    Arguments:
    - left: the left-order of the contraction
    - right: the right-order of the contraction
    - divide_freqs: returns the coefficient as a function of only one frequency array ω

    Returns:
    - c(μ, ν): a symbolic expression for the contraction coeffeicient.
    - c(ω): a symbolic expression for the contraction coeffeicient, with just one argument.
"""
function contraction_coeff(left::Int, right::Int, divide_freqs = false)
    node = DiagramNode((left, right))
    diagrams = get_diagrams(node)
    function _c(μ, ν)
        ex = 0                                                                  # initialize an empty expression
        for diagram in diagrams
            #Check if the diagram has any poles
            s_list, stag_list = find_all_poles(ω)
            total_num_poles = sum([length(su) + length(sl) for (su, sl) in (s_list, stag_list)])
            if total_num_poles == 0 #No singularities means we don't need to Taylor expand
                diag_ex = 1
                for (i, bubble) in enumerate(diagram)
                    endmode = (i == length(diagram)) ? 1 : 0
                    diag_ex = diag_ex*(bubble_coeff(left, right, endmode))(μ, ν)    # calculate the bubble correction term and multiply it in
                end

            else #We calculate the finite term of the Taylor expansion
                @cnumbers τ
                f(x) = exp(-1/2*τ^2*x^2)
                ω = diagram
                # first sum (n+u+l)_i = total_num_poles
                unl_list = find_integer_solutions(3*length(ω), total_num_poles)
                for (i, (μ, ν)) in enumerate(ω)
                    l = length(μ)
                    r = length(ν)
                    
                    su_list = find_poles(μ)
                    sv_list = find_poles(ν)                         # singular indices

                    ju_list = filter(x -> !(x in su_list), 1:l)
                    jl_list = filter(x -> !(x in sv_list), 1:r)     # non-singular indices

                    num_poles = length(su_list) + length(sv_list)   # total number of poles
                
                    # finite part of the bubble factor (not including the expansion terms)
                    A = -f(sum(μ) + sum(ν))/(vec_factorial(μ, include_poles = false)*vec_factorial(ν, include_poles=false))
                
                    unl_list = reshape_sols_vec(find_integer_solutions(3, num_poles), 3, num_poles, length(ω))
                    
                    for (idx, (n, u, d)) in enumerate(unl_list[i, :])
                        numerator_terms = [taylor_coeff(n, k)/factorial(n)*τ^(2*(n - k))*(sum(μ) + sum(ν))^(n - 2*k)*(d + r)^n for k = 1:floor(Int, n/2)]

                        mu_list = find_integer_solutions(length(ju_list), u)
                        ml_list = find_integer_solutions(length(jl_list), d)

                        ###
                        for (mu, ml) in product(mu_list, ml_list)
                            for (ju, jl) in product(ju_list, jl_list)
                                non_poles_product = (-ju/sum(μ(1:ju)))^mu_list(ju)*(-jl/sum(ν(1:jl)))^ml_list(jl)
                            end
                            
                            poles_product = prod[su*sl for (su, sl) in product(su_list, sl_list)]
                            denominator_terms = non_poles_product/poles_product
                        end
                        ###
                        diagram_factors[idx] = A*sum(numerator_terms)*sum(denominator_terms)
                    end
                end
                diag_ex = prod(diagram_factors)
            end
                ex = ex + diag_ex                                                   # add the diagram to the contraction coefficient expression
        end
        return ex
    end
    if divide_freqs
        _c(ω) = (-1)^(left)*_c(ω[1:left], ω[left+1:end])  
    end
    return _c
end


function contraction_coeff(order::Tuple{Int, Int})
    left, right = order
    return contraction_coeff(left, right)
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
