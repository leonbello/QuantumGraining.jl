using Symbolics
using SymbolicUtils
using QuantumGraining

"""
    expand_operators(hs)

Expand a list of operators by applying the `expand_operator` function to each operator.

# Arguments
- `hs`: A list of operators to be expanded.

# Returns
A list of expanded operators.

"""
function expand_operators(hs)
    unique_hs = []

    for h in hs
        new_h, _ = expand_operator(h)
        push!(unique_hs, new_h...)
    end

    return unique_hs
end

"""
    expand_operators(hs, gs, ωs)

Goes over any sum of operators and breaks it into the constituent operators while preserving the order of the other vectors.

# Arguments
- `hs`: Array of operators to be expanded.
- `gs`: Array of coefficients corresponding to the operators.
- `ωs`: Array of weights corresponding to the operators.

# Returns
- `unique_hs`: Array of expanded operators.
- `unique_gs`: Array of coefficients corresponding to the expanded operators.
- `unique_ωs`: Array of weights corresponding to the expanded operators.
"""
function expand_operators(hs, gs, ωs)
    unique_hs = []
    unique_gs = []
    unique_ωs = []

    for (h, g, ω) in zip(hs, gs, ωs)
        new_h, facs = expand_operator(h)

        push!(unique_hs, new_h...)
        push!(unique_gs, [fac*g for fac in facs]...)
        push!(unique_ωs, ω*ones(length(new_h))...)
    end

    return unique_hs, unique_gs, unique_ωs
end

"""
    expand_operator(h)

Expand a quantum operator into a list of individual operators and their corresponding coefficients.

# Arguments
- `h`: The quantum operator to be expanded.

# Returns
- `ops`: A list of individual operators.
- `facs`: A list of corresponding coefficients.

"""
function expand_operator(h)
    ops = []
    facs = []
    if h ≠ 0
        if h isa QuantumCumulants.QAdd  # if it's a sum of operators
            for op in h.arguments
                push!(ops, (op isa QuantumCumulants.QMul) ? prod(op.args_nc) : op)
                push!(facs, (op isa QuantumCumulants.QMul) ? op.arg_c : 1)
            end
        elseif h isa QuantumCumulants.QMul
            push!(ops, h)
            push!(facs, 1)
        end
    end
    return ops, facs
end

"""
    group_operators(hs, gs, ωs; as_dict=true)

Group operators based on their values and return the grouped operators.

# Arguments
- `hs`: An array of operators.
- `gs`: An array of coefficients corresponding to the operators.
- `ωs`: An array of corresponding elements.
- `as_dict`: A boolean indicating whether to return the grouped operators as dictionaries. Default is `true`.

# Returns
- If `as_dict` is `true`, returns two dictionaries `gs_dict` and `ωs_dict` where the keys are the operators and the values are the grouped coefficients and elements respectively.
- If `as_dict` is `false`, returns three arrays `new_hs`, `new_gs`, and `new_ωs` where `new_hs` contains the grouped operators, `new_gs` contains the grouped coefficients, and `new_ωs` contains the grouped elements.

"""
function group_operators(hs, gs, ωs; as_dict=true)
    new_hs = []
    new_gs = []
    new_ωs = []

    for i in eachindex(hs)
        idx = findfirst(x -> isequal(x, hs[i]), new_hs)
        if !isnothing(idx)
            # If it's a duplicate, find its index in new_hs i=14
            new_gs[idx] += gs[i]
        else
            # If it's not a duplicate, add hs[i] to new_hs
            push!(new_hs, hs[i])
            # Add gs[i] to new_gs and add corresponding element from ws to new_ws
            push!(new_gs, gs[i])
            push!(new_ωs, ωs[i])
        end
    end

    new_gs = simplify_contraction.(new_gs)

    if as_dict
        gs_dict = Dict(new_hs[i] => new_gs[i] for i in eachindex(new_hs))
        ωs_dict = Dict(new_hs[i] => new_ωs[i] for i in eachindex(new_hs))
        return gs_dict, ωs_dict
    else
        return new_hs, new_gs, new_ωs
    end
end

"""
    group_operators(hs)

Group operators in the given array `hs` by removing duplicates.

# Arguments
- `hs`: An array of operators.

# Returns
- `new_hs`: An array of operators with duplicates removed.

"""
function group_operators(hs)
    new_hs = []

    for i in eachindex(hs)
        idx = findfirst(x -> isequal(x, hs[i]), new_hs)
        if isnothing(idx)
            push!(new_hs, hs[i])
        end
    end

    return new_hs
end