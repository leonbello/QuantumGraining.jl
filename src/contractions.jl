using IterTools
using QuantumGraining

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
    diagram_correction(ω)
Gives an expression for the effective diagram correction.
    
    # Arguments
    - `ω`: A list of vectors (μi, νi) with the frequency values for each mode. 
    For example: For a diagram d = [(3, 2), (2, 1)] we would have ω = [(μ1, ν1), (μ2, ν2)] where μi and νi are vectors containing mode values.

    # Returns
    - A Correction struct representing the effective diagram correction.
    
 """
 ## New and corrected
function diagram_correction(diagram::Diagram{T1, T2}) where {T1, T2}
    # calculate an array of simple factor terms -- taylor_factors
    simple_factors = calc_simple_factors(diagram)
    total_correction = -prod(simple_factors)

    # if there are poles we need to calculate
    if diagram.num_poles != 0
        num_bubbles = length(diagram)
        sols = find_integer_solutions(3*num_bubbles, diagram.num_poles)    
        sols = reshape_sols(sols, diagram.num_poles, num_bubbles)

        order = 2*diagram.num_poles + 1
        total_poly = zeros(Number, order)
        for i in 1:size(sols)[1]
            poly = zeros(Number, order)
            poly[1] = 1
            total_prefac = zeros(Number, num_bubbles)

            for idx in 1:num_bubbles
                μ, ν = diagram[idx].up, diagram[idx].down
                n = sols[i, idx][1]
                u = sols[i, idx][2]
                d = sols[i, idx][3]
                poly = poly_multiplication(poly, calc_expansion_factors(μ, ν, order, n), order)
                prefac = calc_pole_corrections(μ, ν, μ.poles, ν.poles, u, d)
                total_prefac[idx] = isempty(prefac) ? 0 : sum(prefac)
            end

            total_poly += prod(total_prefac).*poly
            for j in 1:length(total_poly)
                total_poly[j] = simplify(total_poly[j])
            end
        end
        total_correction = extend_correction(total_correction, total_poly)
    end
    return total_correction
end 
function diagram_correction(ω::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
    d = Diagram(ω)
    return diagram_correction(d)
end

function calc_simple_factors(d::Diagram{T1, T2}) where {T1, T2}
    taylor_factors = Vector{Correction}()                                                 # array holding the terms of the outer sum
    for b in d
        μ, ν = b.up, b.down
        l = length(b.down)                                                                # number of down modes in the bubble

        exponent = sum(μ)^2 + sum(ν)^2 + 2*sum(μ)*sum(ν)                                  # exponent of the correction factor
        prefac = (-1)^(l + 1)*1//(vec_factorial(μ)*vec_factorial(ν))                      # prefactor of the correction factor
        poly = Num[1,]                                                                    # polynomial representing the correction factor
        push!(taylor_factors, Correction(prefac, exponent, poly))                         # add the correction factor to the array
    end
    #taylor_factors[1] *= -1                                                              # uncomment this line if the first factor needs to be negated
    return taylor_factors
end

function calc_expansion_factors(mu::BVector{T1}, ν::BVector{T2}, order, n::Int) where {T1, T2}
    if mu.special
        μ = mu[2:end]
    else
        μ = mu
    end

    l, r = length(μ), length(ν)

    #order = 2*(length(μ.poles) + length(ν.poles)) + 1
    poly = zeros(Num, order)
    for k in 0:floor(Int, n/2)
        freq_sum = isequal(sum(mu) + sum(ν), 0) && (n - 2*k) == 0 ? 1 : (sum(mu) + sum(ν))
        l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases
        coeff = taylor_coeff(n, k)//factorial(n)*(freq_sum)^(n - 2*k)*(l_plus_r)^n 
        poly[2*(n-k) + 1] = coeff
    end
    return poly
end

"""
    poly_multiplication(poly1::Vector, poly2::Vector, order::Int)

Multiply two polynomials `poly1` and `poly2` of lengths `len1` and `len2` respectively, and return the resulting polynomial of length `order`.

# Arguments
- `poly1::Vector`: The first polynomial.
- `poly2::Vector`: The second polynomial.
- `order::Int`: The length of the resulting polynomial.

# Returns
- `poly_product::Vector`: The resulting polynomial.

"""
function poly_multiplication(poly1::Vector, poly2::Vector, order::Int)
    len1 = length(poly1)
    len2 = length(poly2)
    poly_product = zeros(Number, order)
    for i in 1:len1
        for j in 1:len2
            if !isequal(abs(poly1[i]*poly2[j]), 0.0) && !isequal(abs(poly1[i]*poly2[j]), -0.0)  # weird comparison issues with Num types
                poly_product[(i-1) + (j-1) + 1] = poly1[i]*poly2[j]
            end
        end
    end
    return poly_product
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

pole_fac(v, regular, j, m) = (-j//sum(v[1:regular]))^m

function calc_pole_corrections(mu::BVector{T1}, ν::BVector{T2}, up_poles::Vector{Int}, down_poles::Vector{Int}, u::Int, d::Int) where {T1, T2}
    if mu.special
        μ = mu[2:end]
    else
        μ = mu
    end
    μ = reverse(μ)

    norm = calc_pole_normalization(up_poles, down_poles)

    pole_terms = []

    up_regular = filter(x -> !(x in up_poles), 1:length(μ))
    down_regular = filter(x -> !(x in down_poles), 1:length(ν))  # non-singular indices
    # change to a more informative name like up_regular and down_regular

    # mu_list and ml_list hold vectors of mu values
    mu_list = find_integer_solutions(length(up_regular), u)
    ml_list = find_integer_solutions(length(down_regular), d)

    if length(up_regular) == 0
        if u == 0
            prod_fac_u = [1]
        else
            prod_fac_u = [0]
        end
    else
        prod_fac_u = []
        for mu_vec in mu_list
            fac_u = []
            for (u, ju) in enumerate(up_regular)
                push!(fac_u, pole_fac(μ, ju, ju, mu_vec[u]))
            end
            push!(prod_fac_u, isempty(fac_u) ? 1 : prod(fac_u))
        end
    end

    if length(down_regular) == 0
        if d == 0
            prod_fac_v = [1]
        else
            prod_fac_v = [0]
        end
    else
        prod_fac_v = []
        for ml_vec in ml_list
            fac_v = []
            for (l, jl) in enumerate(down_regular)
                push!(fac_v, pole_fac(ν, jl, jl, ml_vec[l]))
            end
            push!(prod_fac_v, isempty(fac_v) ? 1 : prod(fac_v))
        end
    end

    for (prod_u, prod_v) in product(prod_fac_u, prod_fac_v)
        push!(pole_terms, prod_u*prod_v//norm)
    end

    return pole_terms   
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
function calculate_bubble_factor(b::Bubble{T1, T2}, sols, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
    μ, ν = b.up, b.down
    exponent = sum(μ)^2 + sum(ν)^2 + 2*sum(μ)*sum(ν) 
    prefac =  1//(vec_factorial(μ)*vec_factorial(ν))
    poly = Float64[1,]
    if !isempty(up_poles) || !isempty(down_poles)
        poly = singular_expansion(b, sols, up_poles, down_poles)
    end
    correction = Correction(prefac, exponent, poly)
    return correction
end

"""
    singular_expansion(b, sols, s, stag; first_bubble = false)
Calculates the taylor expansion coefficients for a singular bubble.

# Arguments
- `B::Bubble{T1, T2}`: the bubble for which we want to calculate the correction factor`
- `sols::Vector{Tuple{Int, Int, Int}}`: vector of solutions to the partition problem
- `s::Vector{Int}`: vector of singular indices in the up-bubbles
- `stag::Vector{Int}`: vector of singular indices in the down-bubbles
- `first_bubble::Bool`: true if this is the first bubble in the chain

# Returns
- `terms::Vector{Complex{Float64}}`: vector of Taylor expansion coefficients
    
"""
function singular_expansion(b::Bubble{T1, T2}, sols, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
    order = 2*(length(up_poles) + length(down_poles)) + 1
    terms = zeros(order)

    for (n, u, d) in sols
        # first inner sum -- includes all analytic contributions
        taylor_terms = calc_analytic_terms(b, order, n)      

        # second inner sum -- factor due to finite pole contributions for this particular solution
        correction_terms = calc_pole_factor(b, u, d, up_poles, down_poles)    
        correction_sum = isempty(correction_terms) ? 0 : sum(correction_terms)

        terms = [terms[i] + correction_sum*taylor_terms[i] for i in 1:order]
    end
    return terms
end

"""
    calc_analytic_terms(b::Bubble{T1, T2}, n)
Calculates the analytic part of the expansion for the diagram, returns an array representing a polynomial.

# Arguments
- `b::Bubble{T1, T2}`: vector of up mode frequencies
- `n::Int`: index for the sum
"""
function calc_analytic_terms(b::Bubble{T1, T2}, order, n) where {T1, T2}
    μ, ν = b.up, b.down
    l, r = length(μ), length(ν)

    # freq_sum is not over all frequencies, just those that do not sum up to 0
    poly = zeros(Float64, order)
    for k in 0:floor(Int, n/2)
        freq_sum = isequal(sum(μ) + sum(ν), 0) && (n - 2*k) == 0 ? 1 : (sum(μ) + sum(ν))
        l_plus_r = (l + r) == 0 && n == 0 ? 1 : (l + r)     # explicity deals with the 0^0 cases
        coeff = taylor_coeff(n, k)//factorial(n)*(freq_sum)^(n - 2*k)*(l_plus_r)^n 

        poly[2*(n-k) + 1] += coeff
    end
    return poly
end

function calc_pole_factor(b::Bubble{T1, T2}, u::Int, d::Int, up_poles::Vector{Int}, down_poles::Vector{Int}) where {T1, T2}
    pole_terms = []
    μ, ν = b.up, b.down
    l, r = length(μ), length(ν)

    ju_list = filter(x -> !(x in up_poles), 1:l)
    jl_list = filter(x -> !(x in down_poles), 1:r)  # non-singular indices

    # mu_list and ml_list hold vectors of mu values
    mu_list = find_integer_solutions(length(ju_list), u)
    ml_list = find_integer_solutions(length(jl_list), d)

    # denominator normalization factor - equals to 1 if s or s' is empty.
    denominator = begin
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

    for (mu_vec, ml_vec) in product(mu_list, ml_list)   
        fac_u = []
        for ju in ju_list
            idx_u = indexin(ju, ju_list)[1]                 # pick up the corresponding index for m_{Ju}
            push!(fac_u, norm_fac(μ, ju, mu_vec[idx_u]))
        end
        prod_fac_u = isempty(fac_u) ? 1 : prod(fac_u)

        fac_v = []
        for jl in jl_list
            idx_l = indexin(jl, jl_list)[1]                 # pick up the corresponding index for m_{Jl}
            push!(fac_v, norm_fac(ν, jl, ml_vec[idx_l]))
        end
        prod_fac_v = isempty(fac_v) ? 1 : prod(fac_v)

        push!(pole_terms, prod_fac_u*prod_fac_v//denominator) 
    end
    return pole_terms
end
