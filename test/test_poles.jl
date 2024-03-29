using Symbolics
using SymbolicUtils
using IterTools

module Tst
    using Test
    using IterTools
    include("../src/bvector.jl")
    include("../src/bubble.jl")
    include("../src/diagram.jl")
    include("../src/diagrams.jl")
    include("../src/poles.jl")
    include("../src/printing.jl")

    #===
    find_poles()
    ===#
    #@cnumbers a b c

    # Using QuantumCumulants for symbolic variables is generally not a good idea, 
    # since it creates type inconsistencies between the vectors. Should only switch right before having to use the operator algebra capabilities.
    @variables a b c
    μ1 = [a, b, c]
    ν1 = [2*a, -2*a]
    ω = [(μ1, ν1)]

    @show find_poles(μ1[2:end])
    @show find_poles(ν1)
    
    #===
    find_all_poles()
    ===#
    μ1 = [0, 1, -1] # pole at 3 -- should be regarded as 2 since the first mode is omitted
    ν1 = [1, -2]
    ω = [(μ1, ν1)]

    find_poles(μ1[2:end])
    
    @show typeof(ω)
    s_list, stag_list = find_all_poles(ω)
    count_poles(s_list, stag_list)
    
    # non-singular indices
    s = s_list[1]
    stag = stag_list[1]
    ju_list = filter(x -> !(x in s), 1:3)      # should be [1, 3]
    jl_list = filter(x -> !(x in stag), 1:2)   # should be [1, 2]
    
    # non-singular cases
    μ1 = [a, b, -a-b] # pole at 3
    ν1 = [a, -a] # pole at 2

    μ2 = [a, a, -2*a] # pole at 1 and 4

    ν2 = [-a, b, a-b] # poles at 3
    
    (μ1, ν1)
    (μ2, ν2)
    ω = [(μ1, ν1), (μ2, ν2)]

    # problem with bubble vectors being different types... need to find a robust way of handling this
    s_list, stag_list = find_all_poles(ω)
    total_num_poles = sum([length(su) + length(sl) for (su, sl) in (s_list, stag_list)])

    @show s_list
    @show stag_list
    @show total_num_poles

    ## taylor_coeff() ##
    # base cases verification
    @show taylor_coeff(1, -1)
    @show taylor_coeff(231, -1)
    @show taylor_coeff(3, 5)
    @show taylor_coeff(3, 3)
    println("---------")

    # general cases
    for i in 1:6
        for j in 1:6
            val = taylor_coeff(i, j)
            if val != 0  
                println("taylor_coeff($i, $j) = $(val)")
                println()
            end
        end
    end

    # comparison to Wentao's code -- works
 
    #=
        find_integer_solutions() 
    weird mismatch between the number of solutions I expect and the number of solutions I get
    =#
    @show find_integer_solutions(1, 5)
    @show find_integer_solutions(2, 4)
    num_indices = 3
    num_bubbles = 2
    num_vars = num_indices*num_bubbles
    target_sum = 6

    num_sols = binomial(target_sum + num_vars - 1, num_vars - 1)
    sols = find_integer_solutions(num_vars, target_sum)
    

    #= 
        reshape_sols() 
    =#
    vecs = reshape_sols(sols, target_sum, num_bubbles)


    @show size(vecs) 
    @show num_sols
    @show vecs

    
end