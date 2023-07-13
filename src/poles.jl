using QuantumCumulants
using SymbolicUtils

"""
    taylor_coeff(n, k)

Recursive function used to calculate finite contribution for Taylor expansion for diagrams with poles.

"""
function taylor_coeff(n::Int, k::Int)
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
find_poles(u)

Finds all the factors of the vector factorial that evaluate to 0.

Argument: 
    - ω = [(μ_1, ν_1), ..., (μ_||d||, ν_||d||)]

Returns:
    - poles_list: list of indices of poles
"""
function find_poles(u::Vector{T}) where {T}
    poles_list = Int[]
    for i = 1:length(u)
        temp_sum = sum(u[1:i])
        if isequal(temp_sum, 0)
            append!(poles_list, i)
        end
    end
    return poles_list
end
find_poles(u::BVector) = u.special ? find_poles(u.freqs[2:end]) : find_poles(u.freqs)



"""
    find_all_poles(d::Vector{Tuple{BVector, BVector}})
Finds all vector factorial poles for a list of BVectors. 
"""
function find_all_poles(d::Vector{Tuple{BVector{T1}, BVector{T2}}}) where {T1, T2}
    up_poles = Vector{Vector{Int}}()
    down_poles = Vector{Vector{Int}}()
    for (μ, ν) in d
        push!(up_poles, find_poles(μ))
        push!(down_poles, find_poles(ν))
    end
    return (up_poles, down_poles)
end

"""
    find_all_poles(d::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
Finds all poles in vector factorial for frequency list `d`. Assumes that the first tuple in the list is the special mode.

Argument: 
    - ω = [(μ_1, ν_1), ..., (μ_||d||, ν_||d||)]
Returns:
    - up_poles: list of indices of poles in upper modes
    - down_poles: list of indices of poles in lower modes
"""
function find_all_poles(freqs::Vector)#::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
    d = Diagram(freqs)
    return find_all_poles(d)
end
find_all_poles(d::Diagram) = find_all_poles(d.freqs)




"""
    count_poles(s_list::Vector{Int}, stag_list::Vector{Int})

Given two lists of poles, counts the total number of poles by counting the non-empty lists.

Argument: 
    - s_list: list of upper poles
    - stag_list: list of lower poles
"""
function count_poles(s_list::Vector{Vector{Int}}, stag_list::Vector{Vector{Int}})
    count = 0
    for list in (s_list, stag_list)
        if !isempty(list)
            for vec in list
                count += length(vec)
            end
        end
    end
    return count
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
function reshape_sols(sols, target_sum, num_bubbles, num_indices=3)
    num_vars = num_bubbles*num_indices
    num_sols = binomial(target_sum + num_vars - 1, num_vars - 1)    

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


# function find_all_poles(freqs::Vector{Tuple{Vector{T1}, Vector{T2}}}) where {T1, T2}
#     up_poles = Vector{Vector{T1}}()
#     down_poles = Vector{Vector{T2}}()
#     for (idx, (μ, ν)) in enumerate(freqs)
#         start = (idx == 1) ? 2 : 1                  # omit the first mode in the first bubble
#         push!(up_poles, find_poles(μ[end:-1:start]))
#         push!(down_poles, find_poles(ν[start:end]))
#     end
#     return (up_poles, down_poles)
# end