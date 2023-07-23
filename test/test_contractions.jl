using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
#using Test

#@testset "contractions" begin
module Tst
    using Test
    using IterTools
    include("../src/diagrams.jl")
    include("../src/bvector.jl")
    include("../src/bubble.jl")
    include("../src/diagram.jl")
    include("../src/poles.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")
    include("../src/corrections.jl")

    ## problem definition ##    
    num_bubbles = 2
    @cnumbers a b c
    μ1 = [0, a, -a] # pole at 1 and 3
    ν1 = [a, -a, b, c] # pole at 2

    μ2 = [a, b, c, -(a+b+c)] # pole at 4
    ν2 = [3*c, 2*b] # no poles

    ω = [(μ1, ν1), (μ2, ν2)]
    d = Diagram(ω)

    s_list = [[1, 3], [4]]
    stag_list = [[2], []]  # singular indices
    total_num_poles = 4
    
    num_vars = 3*num_bubbles
    num_sols = binomial(total_num_poles + num_vars - 1, num_vars - 1)
    sols = find_integer_solutions(3*num_bubbles, total_num_poles)     
    unl_list = reshape_sols(sols, total_num_poles, num_bubbles)           # partition for the inner sum

    ## singular_expansion & calculate_bubble_factor() - bubble 1 ##
    s = s_list[1]
    stag = stag_list[1]
    
    @show s, stag
    t1 = singular_expansion(d[1], unl_list[:, 1], s, stag)
    t2 = calculate_bubble_factor(d, 1, unl_list[:, 1], s, stag)

    #Testing ContractionCoefficient 
    
    μ1 = [1, 1, 3];  
    ν1 = [1, -2];
    μ2 = [3, 1];
    ν2 = Int[];
    ω = Diagram([(μ1, ν1), (μ2, ν2)]);

    μ1 = [0, 1, -2];  
    ν1 = [];
    μ2 = [1, 3];
    ν2 = [];
    ω2 = Diagram([(μ1, ν1), (μ2, ν2)]);

    μ1 = [0, 2, -2];  
    ν1 = [];
    μ2 = [1, 3];
    ν2 = [];
    ω3 = Diagram([(μ1, ν1), (μ2, ν2)]);

    corr1 = diagram_correction(ω)
    corr2 = diagram_correction(ω2)
    corr3 = diagram_correction(ω3)
    typeof(corr1.poly)
    c1 = corr1 + corr2 
    c1.prefacs
    c2 = c1 + corr2
    contract2.exponents
    contract2.prefacs
    contract = c1 + contract2
    contract.exponents

    corr1 + corr2 + corr2 + corr1 + corr3

end



#=
end #testset
=#