using IterTools

"""
    split_freqs_into_bubbles(freqs, diagram)
Splits an array of frequencies into an array of tuples of frequencies, matching the dimensions of each bubble in `diagram`.

Argument:
    freqs - an array of frequencies in the form ω = [μ1..., μ2..., ..., μl..., ν1..., ν2..., ..., νr]
    diagram - dimensions of each bubble, where each tuple is the dimension for the corresponding bubble.
Returns:
    ω - an array of frequencies in the form [(μ1, ν1), (μ2, ν2), ..., (μl, νl), ([], ν(l+1)), ..., ([], νr)]
"""
function split_freqs_into_bubbles(freqs::Vector, diagram::Vector{Tuple{Int, Int}})
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
    contraction_coeff(left::Int, right::Int, freqs::Array)

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
        reverse!(diagram)                               # reversing since Wentao's order is right-to-left, rather than left-to-right
        ω = split_freqs_into_bubbles(freqs, diagram)
        push!(c_list, (diagram_correction(ω), diagram))
        c += c_list[end][1] 
    end
    return c, c_list
end
contraction_coeff(order::Tuple{Int, Int}, ω::Array) = contraction_coeff(order[1], order[2], ω)


 ## New and corrected
function diagram_correction(diagram::Diagram{T1, T2}) where {T1, T2}
    """
        diagram_correction(diagram::Diagram)
    Gives an expression for the effective diagram correction.
    
    # Arguments
    - `diagram`: A diagram object containing all the mode frequencies.

    # Returns
    - A Correction struct representing the diagram correction.
    
    """
    # calculate an array of simple factor terms 
    simple_factors = calc_simple_factors(diagram)
    total_correction = prod(simple_factors)

    if diagram.num_poles != 0                                                   # if there are poles we need to calculate
        num_bubbles = length(diagram)
        sols = find_integer_solutions(3*num_bubbles, diagram.num_poles)    
        sols = reshape_sols(sols, diagram.num_poles, num_bubbles)
        
        ## debugging: make sure the partitions are correct.
        # @show num_bubbles
        # @show sols
        ##

        order = 2*diagram.num_poles + 1
        total_poly = zeros(Number, order)
        for (i, (n, u, d)) in enumerate(sols)
            bubble_idx = (i - 1) % num_bubbles + 1                              # repeats: num_bubbles + i → i
            
            ## debugging: make sure the partitions are correct
            # @show i, bubble_idx
            # @show (n, u, d)
            ##

            μ, ν = diagram[bubble_idx].up, diagram[bubble_idx].down
            poly = calc_expansion_factors(μ, ν, order, n)                       # calculates a list of polynomial coeffecients for a given `n`
            prefac = calc_pole_corrections(μ, ν, μ.poles, ν.poles, u, d)
            if !isempty(prefac)
                total_poly += prefac.*poly
            end
        end
        total_correction = extend_correction(total_correction, total_poly)
    end
    return total_correction
end

function diagram_correction(ω::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
    """
        diagram_correction(ω::Vector{Tuple{Vector{T1}, Vector{T2}}})
    An onverloaded version of `diagram_correction` for vectors of frequencies..

    # Arguments
    - `ω`: A list of vectors (μi, νi) with the frequency values for each mode. 
    For example: For a diagram d = [(3, 2), (2, 1)] we would have ω = [(μ1, ν1), (μ2, ν2)],
    where μi and νi are vectors containing mode values.
    """
    d = Diagram(ω)
    return diagram_correction(d)
end

function calc_simple_factors(diagram::Diagram{T1, T2}) where {T1, T2}
    """
        calc_simple_factors(d::Diagram)
    Calculates the simple analytic factors for a diagram contribution.
    ``
        \\mathcal{C}(\\mu, \\nu) = (-1)^l \\prod_{i=1}^{||d||} \\frac{f(\\sum_i \\mu_i + \\nu_i)}{\\vec{\\mu}!\\vec{\\nu}!}
    ``

    # Arguments
    - `d`: A diagram object including all the different frequency modes for each bubble μᵢ and νᵢ.

    # Returns
    - A Correction struct representing the diagram correction due to only the analytical terms.
    """
    factors = Vector{Correction}()                                       # array holding the terms of the outer sum
    for bubble in diagram
        μ, ν = bubble.up, bubble.down
        l = length(bubble.up)

        exponent = sum(μ)^2 + sum(ν)^2 + 2*sum(μ)*sum(ν) 
        prefac =  (-1)^l*1/(vec_factorial(μ)*vec_factorial(ν))
        poly = Num[1,]
        push!(factors, Correction(prefac, exponent, poly))
    end
    return factors
end


function calc_expansion_factors(μ::BVector{T1}, ν::BVector{T2}, order, n::Int) where {T1, T2}
    """
        calc_expansion_factors(μ::BVector{T1}, ν::BVector{T2}, order, n::Int) where {T1, T2}
    Calculates the terms of the second inner sum in an array holding the coefficients of the power of τ.
    ``
        P[c] = \\sum_{2(n_i - k_i)=c}^{\\floor{n_i/2}} \\frac{c(n_i, k_i)}{n_i!} (\\sum_i (\\mu_i + \\nu_i))^{n_i - 2 k_i} (l + r)^{n_i}
    ``
        
        # Arguments
        - `μ`, `ν``: vectors of the modes for this bubble.
        -  `order``: the order of the polynomial (the maximum power of τ).
        - `n`: current `n` index for the outer sum.

        # Returns
        - A list of polynomial coefficients.
        
    """
    l, r = length(μ), length(ν)
    poly = zeros(Num, order)
    for k in 0:floor(Int, n/2)
        freq_sum = isequal(sum(μ) + sum(ν), 0) && (n - 2*k) == 0 ? 1 : (sum(μ) + sum(ν))
        l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases

        poly[2*(n-k) + 1] = taylor_coeff(n, k)/Float64(factorial(n))*(freq_sum)^(n - 2*k)*(l_plus_r)^n 
    end
    return poly
end
    
function calc_pole_normalization(up_poles, down_poles)
    if isempty(up_poles) && isempty(down_poles)
        return 1
    elseif isempty(up_poles) && !isempty(down_poles)
        return prod(down_poles)
    elseif isempty(down_poles) && !isempty(up_poles)
        return prod(up_poles)
    elseif !isempty(up_poles) && !isempty(down_poles)
        return prod(up_poles)*prod(down_poles)
    end
end

function create_dict(m::Vector{T}, j::Vector{Int}) where T
    dict = Dict{Int, T}()
    for (idx, ji) in enumerate(j)
        dict[ji] = m[idx]
    end
    dict
end

pole_fac(v, j, m) = (-j/sum(v[1:j]))^m


function calc_pole_corrections(μ::BVector{T1}, ν::BVector{T2}, 
    up_poles::Vector{Int}, down_poles::Vector{Int}, 
    u::Int, d::Int) where {T1, T2}
    """
        calc_pole_corrections(μ::BVector{T1}, ν::BVector{T2}, 
                            up_poles::Vector{Int}, down_poles::Vector{Int}, 
                            u::Int, d::Int)
    Calculates the last sum, which scales the polynomial coefficients for each bubble.
    ``
    \\matchal{E}_{u,l}(\\mu, \\nu) = \\frac{1}{\\mathcal{N}} \\prod_{J_u ∉ s_i, J_l ∉ s_i'} 
    \\left \\frac{-J_u}{\\mu^1 + \\cdots \\mu^{J_u}} \\right
    \\left \\frac{-J_l}{\\nu^1 + \\cdots \\nu^{J_l}} \\right
    ``

    The notation here follows that of the equations in the paper, to 
    - (u, d): partition probelm values
    - ju: list of the regular indices of μ, ju = [ju[i] ∉ up_poles] 
    - jl: list of the regular indices of ν, jl = [jl[i] ∉ down_poles]
    - mu: solutions of the partition Σm_ju = u, divided into vectors
    - ml: solutions of the partition Σm_jl = u, divided into vectors
    """

    norm = calc_pole_normalization(up_poles, down_poles)
    l, r = length(μ), length(ν)

    ju_list = filter(x -> !(x in up_poles), 1:l)
    jl_list = filter(x -> !(x in down_poles), 1:r)  # non-singular indices

    mu_list = find_integer_solutions(length(ju_list), u)
    ml_list = find_integer_solutions(length(jl_list), d)

    ## debugging:
    @show ju_list, mu_list
    @show jl_list, ml_list
    ##
    
    pole_terms = []
    for (mu, ml) in product(mu_list, ml_list)       # for each solution of the partitions
        mu = create_dict(mu, ju_list)               # creates a dictionary such that mu[ju[i]] = mu[i]
        ml = create_dict(ml, jl_list)

        ## debugging:
        # @show mu, ju_list
        # @show ml, jl_list
        ##

        fac_u = [pole_fac(μ, ju, mu[ju]) for ju in ju_list]
        prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

        fac_v = [pole_fac(ν, jl, ml[jl]) for jl in jl_list]
        prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

        push!(pole_terms, prod_fac_u*prod_fac_v/norm) 
    end

    sum_pole_terms = !isempty(pole_terms) ? sum(pole_terms) : 0
    return sum_pole_terms   
end



######################################
### May also need to be deprecated ###
######################################
# function diagram_correction(d::Diagram{T1, T2}) where {T1, T2}
#     num_bubbles = length(d)
#     sols = find_integer_solutions(3*num_bubbles, d.num_poles)    
#     sols = reshape_sols(sols, d.num_poles, num_bubbles)         
#     bubble_factors = Vector{Correction}()                                                 # array holding the terms of the outer sum
    
#     l_tot, r_tot = 0, 0
#     for (i, b) in enumerate(d)
#         fac = calculate_bubble_factor(b, sols[:, i], d.up_poles[i], d.down_poles[i])
#         push!(bubble_factors, fac)
#         l_tot += length(b.up)
#         r_tot += length(b.down)
#     end
#     return (-1)^l_tot*prod(bubble_factors)
# end


# """
#     calculate_bubble_factor(ω, bubble_idx, total_num_poles, s, stag)
# Returns the bubble factor for a single bubble.

# # Arguments
#     - `ω`: the frequency values for the whole diagram.
#     - `total_num_poles`: the total number of singular poles in the diagram
#     - `bubble_idx`: the index of the bubble for which we want to calculate the correction factor.
#     - `s`: a list of the singular poles in the upper modes of the bubble.
#     - `stag`: a list of the singular poles in the lower modes of the bubble.

# # Returns 
#     - an array of all bubble factors.
# """
# function calculate_bubble_factor(b::Bubble{T1, T2}, sols, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
#     μ, ν = b.up, b.down
#     exponent = sum(μ)^2 + sum(ν)^2 + 2*sum(μ)*sum(ν) 
#     prefac =  1/(vec_factorial(μ)*vec_factorial(ν))
#     poly = Float64[1,]
#     if !isempty(up_poles) || !isempty(down_poles)
#         poly = singular_expansion(b, sols, up_poles, down_poles)
#     end
#     correction = Correction(prefac, exponent, poly)
#     return correction
# end

# """
#     singular_expansion(b, sols, s, stag; first_bubble = false)
# Calculates the taylor expansion coefficients for a singular bubble.

# # Arguments
# - `B::Bubble{T1, T2}`: the bubble for which we want to calculate the correction factor`
# - `sols::Vector{Tuple{Int, Int, Int}}`: vector of solutions to the partition problem
# - `s::Vector{Int}`: vector of singular indices in the up-bubbles
# - `stag::Vector{Int}`: vector of singular indices in the down-bubbles
# - `first_bubble::Bool`: true if this is the first bubble in the chain

# # Returns
# - `terms::Vector{Complex{Float64}}`: vector of Taylor expansion coefficients
    
# """
# function singular_expansion(b::Bubble{T1, T2}, sols, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
#     # not sure the order is correct here
#     order = 2*(length(up_poles) + length(down_poles)) + 1
#     terms = zeros(order)

#     for (n, u, d) in sols
#         # first inner sum -- includes all analytic contributions
#         taylor_terms = calc_analytic_terms(b, order, n)      

#         # second inner sum -- factor due to finite pole contributions for this particular solution
#         correction_terms = calc_pole_factor(b, u, d, up_poles, down_poles)    
#         correction_sum = isempty(correction_terms) ? 0 : sum(correction_terms)

#         terms = [terms[i] + correction_sum*taylor_terms[i] for i in 1:order]
#     end
#     return terms
# end

# """
#     calc_analytic_terms(b::Bubble{T1, T2}, n)
# Calculates the analytic part of the expansion for the diagram, returns an array representing a polynomial.

# # Arguments
# - `b::Bubble{T1, T2}`: vector of up mode frequencies
# - `n::Int`: index for the sum
# """
# function calc_analytic_terms(b::Bubble{T1, T2}, order, n) where {T1, T2}
#     μ, ν = b.up, b.down
#     l, r = length(μ), length(ν)

#     #poly = convert(Array{Any}, [1, zeros(2*n)...])
#     # freq_sum is not over all frequencies, just those that do not sum up to 0
#     poly = zeros(Float64, order)
#     for k in 0:floor(Int, n/2)
#         freq_sum = isequal(sum(μ) + sum(ν), 0) && (n - 2*k) == 0 ? 1 : (sum(μ) + sum(ν))
#         l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases
#         coeff = taylor_coeff(n, k)/Float64(factorial(n))*(freq_sum)^(n - 2*k)*(l_plus_r)^n 

#         poly[2*(n-k) + 1] += coeff
#     end
#     return poly
# end

# function calc_pole_factor(b::Bubble{T1, T2}, u::Int, d::Int, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
#     pole_terms = []
#     μ, ν = b.up, b.down
#     l, r = length(μ), length(ν)

#     ju_list = filter(x -> !(x in up_poles), 1:l)
#     jl_list = filter(x -> !(x in down_poles), 1:r)  # non-singular indices

#     # mu_list and ml_list hold vectors of mu values
#     mu_list = find_integer_solutions(length(ju_list), u)
#     ml_list = find_integer_solutions(length(jl_list), d)

#     # denominator normalization factor - equals to 1 if s or s' is empty.
#     denominator = begin
#         if isempty(up_poles) && isempty(down_poles)
#             return 1
#         elseif isempty(up_poles) && !isempty(down_poles)
#             return prod(down_poles)
#         elseif isempty(down_poles) && !isempty(up_poles)
#             return prod(up_poles)
#         elseif !isempty(up_poles) && !isempty(down_poles)
#             return prod(up_poles)*prod(down_poles)
#         end
#     end

#     for (mu_vec, ml_vec) in product(mu_list, ml_list)   
#         fac_u = []
#         for ju in ju_list
#             idx_u = indexin(ju, ju_list)[1]                 # pick up the corresponding index for m_{Ju}
#             push!(fac_u, norm_fac(μ, ju, mu_vec[idx_u]))
#         end
#         prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

#         fac_v = []
#         for jl in jl_list
#             idx_l = indexin(jl, jl_list)[1]                 # pick up the corresponding index for m_{Jl}
#             push!(fac_v, norm_fac(ν, jl, ml_vec[idx_l]))
#         end
#         prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

#         push!(pole_terms, prod_fac_u*prod_fac_v/denominator) 
#     end
#     return pole_terms
# end



# #############################
# ### CODE TO BE DEPRECATED ###
# #############################

# function calc_pole_factor(b::Bubble{T1, T2}, s, stag, mu_list, ml_list, ju_list, jl_list) where {T1, T2}
#     pole_terms = []
#     μ, ν = b.up, b.down
#     l, r = length(μ), length(ν)

#     # denominator normalization factor - equals to 1 if s or s' is empty.
#     denominator = begin
#         if isempty(s) && isempty(stag)
#             return 1
#         elseif isempty(s) && !isempty(stag)
#             return prod(stag)
#         elseif isempty(stag) && !isempty(s)
#             return prod(s)
#         elseif !isempty(s) && !isempty(stag)
#             return prod(s)*prod(stag)
#         end
#     end

#     for (mu_vec, ml_vec) in product(mu_list, ml_list)   
#         fac_u = []
#         for ju in ju_list
#             idx_u = indexin(ju, ju_list)[1]                 # pick up the corresponding index for m_{Ju}
#             push!(fac_u, norm_fac(μ, ju, mu_vec[idx_u]))
#         end
#         prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

#         fac_v = []
#         for jl in jl_list
#             idx_l = indexin(jl, jl_list)[1]                 # pick up the corresponding index for m_{Jl}
#             push!(fac_v, norm_fac(ν, jl, ml_vec[idx_l]))
#         end
#         prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

#         push!(pole_terms, prod_fac_u*prod_fac_v/denominator) 
#     end
#     return pole_terms
# end

# # old version, may need to be deprecated
# function calculate_bubble_factor(ω::Vector{Tuple{Vector{T1}, Vector{T2}}}, bubble_idx::Int, sols, s, stag) where {T1, T2}
#     μ, ν = ω[bubble_idx]   
#     ν = reverse(ν)
#     l = length(μ)
#     r = length(ν)
    
#     @cnumbers τ
#     f(x) = exp(-0.5*τ^2*x^2)
    
#     # finite part of the bubble factor (not including the expansion terms)
#     first_bubble = (bubble_idx == 1)
#     μ0, ν0 = set_indices(l, r, first_bubble)

#     sum_μ = isempty(μ) ? 0 : sum(μ)
#     sum_ν = isempty(ν) ? 0 : sum(ν)  # explicity deal with the case where one the vectors is empty (up-bubble or down-bubble)
#     prefac = -f(sum_μ + sum_ν)/(vec_factorial(μ[end:-1:μ0], include_poles = false)*vec_factorial(ν[end:-1:ν0], include_poles=false))

#     return prefac*sum( singular_expansion( μ[μ0:end], ν[ν0:end], sols, s, stag, first_bubble = first_bubble) )
# end

# """
#     singular_expansion(μ, ν, sols, s, stag; first_bubble = false)
# Calculates the taylor expansion coefficients for a singular bubble.

# # Arguments
# - `μ::Vector{Int}`: vector of up mode frequencies
# - `ν::Vector{Int}`: vector of down mode frequencies
# - `sols::Vector{Tuple{Int, Int, Int}}`: vector of solutions to the partition problem
# - `s::Vector{Int}`: vector of singular indices in the up-bubbles
# - `stag::Vector{Int}`: vector of singular indices in the down-bubbles
# - `first_bubble::Bool`: true if this is the first bubble in the chain

# # Returns
# - `terms::Vector{Complex{Float64}}`: vector of taylor expansion coefficients
    
# """
# function singular_expansion(μ, ν, sols, s, stag; first_bubble = false)
#     l = length(μ)
#     r = length(ν)

#     ju_list = filter(x -> !(x in s), 1:l)
#     jl_list = filter(x -> !(x in stag), 1:r)  # non-singular indices

#     order = length(s) + length(stag)
#     terms = zeros(order)

#     # In each loop we calculate a list of polynomial coefficients, for a given solution of the partition problem
#     # We want to have an array that holds the polynomial coefficients of the sum of all solutions
#     for (idx, (n, u, d)) in enumerate(sols)
#         # first inner sum -- includes all analytic contributions
#         # analytic terms should return a list of terms, each term being a polynomial
#         taylor_terms = calc_analytic_terms(μ, ν, n)      
        
#         # mu_list and ml_list hold vectors of mu values
#         mu_list = find_integer_solutions(length(ju_list), u)
#         ml_list = find_integer_solutions(length(jl_list), d)
        
#         # second inner sum -- terms due to finite pole contributions  
#         # this is a correction factor for this particular solution
#         pole_factors = calc_pole_factor(μ, ν, s, stag, mu_list, ml_list, ju_list, jl_list, first_bubble)     
#         poles_sum = isempty(pole_factors) ? 0 : sum(pole_factors)

#         #analytic_sum = isempty(analytic_terms) ? 0 : sum(analytic_terms)

#         terms = [terms[i] + poles_sum*taylor_terms[i] for i in 1:order]
#         #push!(terms, poles_sum.*taylor_terms)
#     end
#     return terms
# end

# """
#     calc_analytic_terms(μ, ν, n)
# Calculates the analytic terms for a single bubble.

# # Arguments
# - `μ::Vector{Int}`: vector of up mode frequencies
# - `ν::Vector{Int}`: vector of down mode frequencies
# - `n::Int`: index for the sum
# """
# function calc_analytic_terms(μ, ν, n)
#     l, r = length(μ), length(ν)
#     @cnumbers τ
    
#     poly = convert(Array{Any}, [1, zeros(2*n)...])
#     #analytic_terms = []
#     for k in 0:floor(Int, n/2)
#         sum_μ = isempty(μ) ? 0 : sum(μ)  
#         sum_ν = isempty(ν) ? 0 : sum(ν)
        
#         freq_sum = isequal(sum_μ + sum_ν, 0) && (n - 2*k) == 0 ? 1 : (sum_μ + sum_ν)
#         l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases
#         coeff = taylor_coeff(n, k)/Float64(factorial(n))*(freq_sum)^(n - 2*k)*(l_plus_r)^n 

#         poly[2*(n-k) + 1] += coeff

#         #push!(analytic_terms, taylor_coeff(n, k)/Float64(factorial(n))*τ^(2*(n - k))*(freq_sum)^(n - 2*k)*(l_plus_r)^n)
#     end
#     #return analytic_terms
#     return poly
# end