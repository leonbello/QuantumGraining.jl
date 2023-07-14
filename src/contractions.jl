using IterTools

"""
    split_freqs_into_bubbles(freqs, diagram)
Splits an array of frequencies into an array of tuples of frequencies, matching the dimensions of each bubble in `diagram`.

Arguments:
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

# Arguments:
    - left: the left-order of the contraction
    - right: the right-order of the contraction
    - freqs: array of frequencies to put in each mode

# Returns:
    - contraction: a contraction struct that includes a list of all the different diagram corrections in an ordered form
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
    return dict
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