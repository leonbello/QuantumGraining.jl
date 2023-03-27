using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
using Test

#@testset "corrections" begin
module Tst
    using Test
    using IterTools
    include("../src/diagrams.jl")
    include("../src/poles.jl")
    include("../src/contractions.jl")
    include("../src/printing.jl")
    include("../src/expressions.jl")
    include("../src/bvector.jl")
    
    # write tests for all basic possible diagram structures
    ## no singularities ##

    # only up-bubbles
    begin   # [(3, 0), (2, 0)], singularity in the first bubble should be omitted
        μ1 = [0, 1, -2];  
        ν1 = [];
        μ2 = [1, 3];
        ν2 = [];
        ω = [(μ1, ν1), (μ2, ν2)];
        
        s, st = find_all_poles(ω);
        @show s
        @show st
        @show count_poles(find_all_poles(ω)...)
        diagram_correction(ω)     
    end

    # only down-bubbles
    begin   # [(0, 3), (0, 2)], singularity in the first bubble should be omitted
        ν1 = [0, 1, -2];   
        μ1 = [];
        ν2 = [1, 3];
        μ1 = [];
        ω = [(μ1, ν1), (μ2, ν2)];

        s, st = find_all_poles(ω);
        @show s
        @show st
        @show count_poles(find_all_poles(ω)...);
        diagram_correction(ω)     
    end

    # general bubbles
    begin   # single bubble [(3, 1)], singularity in the first bubble should be omitted
        μ1 = [0, 1, 2];    
        ν1 = [5];
        ω = [(μ1, ν1)];
        
        @show count_poles(find_all_poles(ω)...)
        diagram_correction(ω)     
    end

    begin   # two bubbles [(3, 2), (2, 1)], singularity in the first bubble should be omitted
        μ1 = [0, 1, 2]
        ν1 = [7, 4]
        μ2 = [1, 3]
        ν2 = [5]
        ω = [(μ1, ν1), (μ2, ν2)];
        
        @show count_poles(find_all_poles(ω)...)
        diagram_correction(ω)     
    end

    ## singularities in the first bubble ##
    # up-singularity in the first up-bubble 
    begin   # [(3, 0), (2, 0)], single singularity in the first bubble
        bubble = [(3,0), (2,0)];
        μ1 = [0, 1, -1];  
        ν1 = [];
        μ2 = [1, 3];
        ν2 = [];
        ω = [(μ1, ν1), (μ2, ν2)];

        @show count_poles(find_all_poles(ω)...)
        diagram_correction(ω)     
    end


    # down-singularity in the first bubble
    begin
        bubble = [(0, 3), (2,0)];    
    end
    

    # up and down singularities in the first bubble
    begin
        bubble = [(3, 2), (1, 2)];    
    end

    ## singularities in the second bubble ##

    # up-singularity in the second bubble 

    # down-singularity in the second bubble 

    # up and down singularities in the second bubble


    # singularities in both bubbles - up-up, up-down, down-up, down-down


    ## singularities in the first (up) bubble ##



    # singularities in the first (down) bubble






    # only up-bubbles
    μ1 = [1, 4]
    ν1 = []
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    @cnumbers ω_1 ω_2
    μ1 = [2]
    ν1 = [1, -1]
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    # only down-bubbles
    ν1 = [1, 4]
    μ1 = []
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    # singular bubbles
    μ1 = [0, 1, -1] # pole at 2
    ν1 = [5]
    ω = [(μ1, ν1)]
    test = diagram_correction(ω)

    @cnumbers ω_1 ω_2

    μ1 = [ω_1, -ω_1]
    ν1 = []
    μ2 = []
    ν2 = [ω_2]
    ω = [(μ1, ν1), (μ2, ν2)] 
    test = diagram_correction(ω)

    μ1 = [ω_1]
    ν1 = []
    μ2 = []
    ν2 = [-ω_2, ω_2]    
    ω = [(μ1, ν1), (μ2, ν2)] #[([ω_1], []), ([], [ω_2, -ω_2])]
    test = diagram_correction(ω)

    # [(1, 2)]
    μ1 = [0] # pole at 2
    ν1 = [1, 1]
    ω1 = [(μ1, ν1)]
    @show diagram_correction(ω1)
    @show c_list[1]

    # [(0, 1), (1, 1)]
    μ1 = [] # pole at 2
    ν1 = [1]
    μ2 = [0]
    ν2 = [1]

    μ1 = [0] # pole at 2
    ν1 = [1]
    μ2 = []
    ν2 = [1]

    ω2 = [(μ1, ν1), (μ2, ν2)]
    diagram_correction(ω2)
    @show c_list[2]
end