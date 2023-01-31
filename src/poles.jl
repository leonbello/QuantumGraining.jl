using QuantumCumulants
using SymbolicUtils

"""
taylor_coeff(n, k)

Recursive function used to calculate finite contribution for Taylor expansion for diagrams with poles.

"""
function taylor_coeff(n, k)
    # base cases
    if (k == -1) || (n < 2*k)
        return 0
    elseif (n == 0) && (k ==0)
        return 1
    else
        return -taylor_coeff(n - 1, k) + (n - 2*k + 1)*taylor_coeff(n - 1, k - 1)
    end
end

"""
vec_factorial(u; include_poles=true)

Calculates vector factorial for a frequency list u. 

"""

function vec_factorial(u; include_poles=true)
    temp_prod = u[1]
    for i = 2:length(u)
        temp_sum = sum(u[1:i])
        if isequal(temp_sum, 0)
            if include_poles == true
                temp_prod = 0
            end
        else
            temp_prod *= temp_sum
        end
    end
    return temp_prod
end

"""
find_poles(u)

Finds all the factors of the vector factorial that evaluate to 0.

Argument: 
    - ω = [(μ_1, ν_1), ..., (μ_||d||, ν_||d||)]
"""

function find_poles(u)
    poles_list = []
    for i = 1:length(u)
        temp_sum = sum(u[1:i])
        if isequal(temp_sum, 0)
            append!(poles_list, i)
        end
    end
    return poles_list
end

"""
find_all_poles(u)

Finds all poles in vector factorial for frequency list u.

Argument: 
    - ω = [(μ_1, ν_1), ..., (μ_||d||, ν_||d||)]
"""

function find_all_poles(ω)
    μ_poles = []
    ν_poles = []
    for (μ, ν) in ω
        push!(μ_poles, find_poles(μ))
        push!(ν_poles, find_poles(ν))
    end
    return (μ_poles, ν_poles)
end


"""
find_integer_solutions(num_vars::Int, target_sum::Int, combination::Vector{Int}=Vector{Int}(), sum_so_far::Int=0)

Function that calculates combinations of k positive integers adding up to m. 

"""
function find_integer_solutions(num_vars::Int, target_sum::Int, combination::Vector{Int}=Vector{Int}(), sum_so_far::Int=0)
    res = []
    if sum_so_far == target_sum
        if length(combination) == num_vars
            push!(res, combination)
            return res
        end
    elseif sum_so_far > target_sum
        return res
    end
    for i in 0:target_sum
        if length(combination) < num_vars
            res = vcat(res, find_integer_solutions(num_vars, target_sum, [combination..., i], sum_so_far+i))
        end
    end
    return res
end


"""
reshape_sols(sols, target_sum, num_bubbles, num_indices = 3)

Helper function that reshapes integer combinations from find_integer_solutions() into vectors
"""
function reshape_sols(sols, target_sum, num_bubbles, num_indices = 3)
    num_vars = num_bubbles*num_indices
    num_sols = floor(Int, factorial(num_vars + target_sum - 1)/(factorial(num_vars)*factorial(target_sum - 1)))
    dim_sols = num_sols                # total number of solutions
    dim_indices = num_indices          # number of indices - u, n, l
    dim_bubbles = num_bubbles          # number of bubbles 
    
    vectors = Array{Array{Int64,1}, 2}(undef, dim_sols, dim_bubbles)

    for i in 1:dim_sols           
        for j in 1:dim_bubbles    
            end_idx = dim_indices*j < length(sols[i][:]) ? dim_indices*j : length(sols[i][:])
            vectors[i, j] = sols[i][(1 + dim_indices*(j - 1)):end_idx]
        end
    end
    return vectors
end

