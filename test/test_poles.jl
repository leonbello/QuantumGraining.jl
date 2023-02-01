#Test file for poles.jl 
using QuantumCumulants

module Tst
    using Test
    include("../src/diagrams.jl")
    include("../src/poles.jl")
    include("../src/printing.jl")

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
    

    """
    ## vec_factorial() ##
    """
    a = []
    vec_factorial(a)
    
    # numeric
    a = [1, 3, 5, 7]
    @show vec_factorial(a) 
    @show 1*(1+3)*(1+3+5)*(1+3+5+7)

    b = [1, 3, 5, 7, -(1+ 3 + 5 + 7)]
    @show vec_factorial(b)
    @show vec_factorial(b, include_poles=false)

    # symbolic
    @cnumbers a b c d

    v = [a, b, c, -(a+b+c), d]
    @show vec_factorial(v, include_poles=false)
    @show vec_factorial(v, include_poles=true)


    """
    find_poles()
    """
    μ1 = [a, b, c]
    ν1 = [2*a, -2*a]
    ω = [(μ1, ν1)]

    @show find_poles(μ1[2:end], start=2)
    @show find_poles(ν1)
    
    """
    find_all_poles()
    """
    μ1 = [0, 1, -1] # pole at 3 -- should be regarded as 2 since the first mode is omitted
    ν1 = [1, -2]
    ω = [(μ1, ν1)]

    find_poles(μ1[2:end], start=2)
    
    s_list, stag_list = find_all_poles(ω)
    count_poles(s_list, stag_list)    
    
    # non-singular cases
    μ1 = [a, b, -a-b] # pole at 3
    ν1 = [a, -a] # pole at 2

    μ2 = [0, a, a, -2a] # pole at 1 and 4
    ν2 = [-a, b, a-b] # poles at 3
    ω = [(μ1, ν1), (μ2, ν2)]

    s_list, stag_list = find_all_poles(ω)
    total_num_poles = sum([length(su) + length(sl) for (su, sl) in (s_list, stag_list)])

    @show s_list
    @show stag_list
    @show total_num_poles

    ## find_integer_solutions() ##

    @show find_integer_solutions(1, 5)
    @show find_integer_solutions(2, 4)

    k = 3*4
    m = 4
    d = 4
    sols = find_integer_solutions(k, m);
    num_sols = floor(Int, factorial(4 + 6 - 1)/(factorial(4)*factorial(5)));

    ## reshape_sols() ##
    # 4 bubbles, each bubble with 3 indices that sum up to 4.
    num_indices = 3
    num_bubbles = 4
    num_vars = num_indices*num_bubbles
    target_sum = 4

    sols = find_integer_solutions(num_vars, target_sum);
    num_sols = floor(Int, factorial(num_vars + target_sum - 1)/(factorial(num_vars)*factorial(target_sum - 1)))

    vecs = reshape_sols(sols, target_sum, num_bubbles)


    @show size(vecs) 
    @show num_sols
    @show vecs

    
end