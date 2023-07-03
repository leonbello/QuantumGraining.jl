using Revise
using IterTools
using Symbolics
using Test
using QuantumGraining

#@testset "corrections" begin
module Tst
    using Test
    using QuantumGraining

    # one common bubble and one up-bubble
    # no singularities
    begin
        μ1 = [1, 1, 3];  
        ν1 = [1, -2];
        μ2 = [3, 1];
        ν2 = Int[];
        ω = [(μ1, ν1), (μ2, ν2)];
        
        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω) 
        @test corr.prefac ≈ 1//48
        @test corr.exponent ≈ 2*16
    end
    ## CORRECT: 1//48*exp(-16*τ^2) ##

    # only up-bubbles
    # no singularities
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
        corr = diagram_correction(ω)

        @test corr.prefac ≈ -1//24
        @test corr.exponent ≈ 2*8.5
    end
    ## CORRECT! -1//24*exp(-8.5τ^2) ##

    # only down-bubbles
    # no singularities
    begin   # [(0, 3), (0, 2)], singularity in the first bubble should be omitted
        μ1 = [];
        ν1 = [0, 1, -2];   
        μ2 = [];
        ν2 = [1, 3];
        
        ω = [(μ1, ν1), (μ2, ν2)];

        s, st = find_all_poles(ω);
        @show s
        @show st
        @show count_poles(find_all_poles(ω)...);
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8.5
    end
    ## Not sure about this one, need to think ##

    # single common bubble
    # no singularities
    begin   # single bubble [(3, 1)], singularity in the first bubble should be omitted
        μ1 = [0, 1, 2];    
        ν1 = [5];
        ω = [(μ1, ν1)];
        
        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)
        
        @test corr.exponent ≈ 2*32
        @test corr.prefac ≈ -1//30
    end
    ## CORRECT! -1//30*exp(-32τ^2) ##

    # two common bubbles
    # no singularities
    begin
        μ1 = [0, 1, 2]
        ν1 = [7, 4]
        μ2 = [1, 3]
        ν2 = [5]
        ω = [(μ1, ν1), (μ2, ν2)];
        
        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)     

        @test corr.prefac ≈ -1//27720
        @test corr.exponent ≈ 277
        @test corr.poly ≈ [1]
    end
    ## CORRECT! -1/27720*exp(-277/2*τ^2) ##

    # only up-bubbles
    begin   # [(3, 0), (2, 0)], single singularity in the first bubble
        μ1 = [0, 1, -1];  
        ν1 = Int[];
        μ2 = [1, 3];
        ν2 = Int[];
        ω = [(μ1, ν1), (μ2, ν2)];

        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        @test corr.poly ≈ [1, 0, -48]
        @test corr.prefac ≈ 1//144
    end

    # single singularity in the first bubble
    begin   # [(3, 0), (2, 0)], single singularity in the first bubble
        μ1 = [0, 1] #, -1];  
        ν1 = Int[];
        μ2 = [-1, 1, 3];
        ν2 = Int[];
        ω = [(μ1, ν1), (μ2, ν2)];

        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)
        corr
        @test corr.exponent ≈ 2*8
        @test corr.poly ≈ [1, 0, -48]
        @test corr.prefac ≈ 1//144
    end
    ## WRONG COEFFICIENTS AND POLY: 1/144*exp(-8*τ^2)*(1-48τ^2) ##

    # up-singularity in the second bubble 
    begin
        μ1 = [1, 3];
        ν1 = [];
        μ2 = [0, 1, -1];  
        ν2 = [];
        
        ω = [(μ1, ν1), (μ2, ν2)];
    
        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        @test corr.poly ≈ [1, 0, -48]
        @test corr.prefac ≈ 1//144
    end


    # up and down singularities in the first bubble
    begin
        μ1 = [0, 1, -1];  
        ν1 = [1, 2, -3];
        μ2 = [1, 3];
        ν2 = Int[]; 
        
        ω = [(μ1, ν1), (μ2, ν2)];
        @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        @test corr.poly ≈ [1, 0, -48]
        @test corr.prefac ≈ 1//144
    end
    ## WRONG COEFFICIENTS


end