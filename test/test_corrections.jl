using Revise
using IterTools
using Symbolics
using Test
#using QuantumGraining

#@testset "corrections" begin
module Tst
    using Test
    include("../src/bvector.jl")
    include("../src/bubble.jl")
    include("../src/diagram.jl")
    include("../src/diagrams.jl")
    include("../src/contractions.jl")
    include("../src/corrections.jl")
    include("../src/lindblad.jl")
    include("../src/poles.jl")
    include("../src/printing.jl")
    include("../src/utils.jl")
    #using QuantumGraining

    # one common bubble and one up-bubble
    # no singularities
    begin
        μ1 = [1, 1, 3];  
        ν1 = [1, -2];
        μ2 = [3, 1];
        ν2 = Int[];
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);
        
        # @show count_poles(find_all_poles(ω)...)
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
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);
        
        # s, st = find_all_poles(ω);
        # @show s
        # @show st
        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.prefac ≈ -1//24
        @test corr.exponent ≈ 2*8.5
    end
    ## CORRECT! -1//24*exp(-8.5τ^2) ##

    # single common bubble
    # no singularities
    begin   # single bubble [(3, 1)], singularity in the first bubble should be omitted
        μ1 = [0, 1, 2];    
        ν1 = [5];
        ω = Diagram([(μ1, ν1)]);
        
        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)
        
        @test corr.exponent ≈ 2*32
        @test corr.prefac ≈ 1//30
    end
    ## CORRECT! 1//30*exp(-32τ^2) ##

    # two common bubbles
    # no singularities
    begin
        μ1 = [0, 1, 2]
        ν1 = [7, 4]
        μ2 = [1, 3]
        ν2 = [5]
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);
        
        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)     

        @test corr.prefac ≈ -1//27720
        @test corr.exponent ≈ 277
    end
    ## CORRECT! -1/27720*exp(-277/2*τ^2) ##

    # only up-bubbles
    begin   # [(3, 0), (2, 0)], single singularity in the first bubble
        μ1 = [0, 1, -1];  
        ν1 = Int[];
        μ2 = [1, 3];
        ν2 = Int[];
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        # @test corr.poly ≈ [1, 0, -48]
        @test isapprox(corr.poly, [1, 0, -48], atol = 1e-13)
        # @test corr.prefac ≈ 1//144
        @test isapprox(corr.prefac, 1//144, atol = 1e-13)
    end
    # there is error in the last one or two digits, probably resulting from limited machine precision

    # single singularity in the first bubble
    begin   # [(3, 0), (2, 0)], single singularity in the first bubble
        μ1 = [0, 1] #, -1];  
        ν1 = Int[];
        μ2 = [-1, 1, 3];
        ν2 = Int[];
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*5
        @test corr.poly ≈ [1]
        @test corr.prefac ≈ -1//36
    end
    ## WRONG COEFFICIENTS AND POLY: 1/144*exp(-8*τ^2)*(1-48τ^2) ##

    # up-singularity in the second bubble 
    begin
        μ1 = [1, 3];
        ν1 = [];
        μ2 = [0, 1, -1];  
        ν2 = [];
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);
    
        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        # @test corr.poly ≈ [1, 0, -69//7, 0, 72//7]
        @test isapprox(corr.poly, [1, 0, -69//7, 0, 72//7], atol = 1e-13)
        # @test corr.prefac ≈ 7//162
        @test isapprox(corr.prefac, 7//162, atol = 1e-13)
    end
    # there is error in the last one or two digits, probably resulting from limited machine precision

    # up and down singularities in the first bubble
    begin
        μ1 = [0, 1, -1];  
        ν1 = [1, 2, -3];
        μ2 = [1, 3];
        ν2 = Int[]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        # @show count_poles(find_all_poles(ω)...)
        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*8
        # @test corr.poly ≈ [1, 0, -90//91, 0, 1152//91]
        @test isapprox(corr.poly, [1, 0, -90//91, 0, 1152//91], atol = 1e-13)
        # @test corr.prefac ≈ 91//7776
        @test isapprox(corr.prefac, 91//7776, atol = 1e-13)
    end



    begin
        μ1 = [0, 1, -1];  
        ν1 = [7, 8];
        μ2 = [1, 3];
        ν2 = [6]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        corr = diagram_correction(ω)

        @test corr.exponent ≈ 325
        @test isapprox(corr.poly, [1, 0, 9450//29], atol = 1e-12)
        @test isapprox(corr.prefac, -29//1587600, atol = 1e-13)
    end

    begin
        μ1 = [1, 2, 3];  
        ν1 = [7, 8];
        μ2 = [4, 5];
        ν2 = [6]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        corr = diagram_correction(ω)

        @test corr.exponent ≈ 2*333
        @test isapprox(corr.poly, [1], atol = 1e-13)
        @test isapprox(corr.prefac, -1//425250, atol = 1e-13)
    end

    begin
        μ1 = [0, 1, -1];  
        ν1 = [0, 0];
        μ2 = [2, -2];
        ν2 = [7]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        corr = diagram_correction(ω)

        @test corr.exponent ≈ 49
        @test isapprox(corr.poly, [1, 0, -2080540//65361, 0, 29868440//65361, 0, -186356016//65361, 0, 311299254//65361], atol = 1e-11)
        @test isapprox(corr.prefac, -65361//4302592, atol = 1e-13)
    end

    begin
        μ1 = [0, 0, 0];  
        ν1 = [2, -2];
        μ2 = [0, 0, 0];
        ν2 = [7]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        corr = diagram_correction(ω)

        @test corr.exponent ≈ 49
        @test isapprox(corr.poly, [1, 0, -58617720//7410735, 0, 38031840//7410735, 0, 38400633600//7410735, 0, 7216608483840//7410735, 0, -107602731651072//7410735, 0, 226775649501184//7410735], atol = 1e-8)
        @test isapprox(corr.prefac, 7410735//113846584320, atol = 1e-13)
    end

    begin
        μ1 = [1, 1, -2, 3, -1, -2, 0];  
        ν1 = Int[];
        μ2 = Int[];
        ν2 = [1, -2, 1, 3, -1, -2, 1, -1]; 
        
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);

        corr = diagram_correction(ω)

        @test corr.exponent ≈ 0
        @test isapprox(corr.poly, [1, 0, -451800//4597019, 0, -1620000//4597019, 0, 0, 0, 0, 0, 0], atol = 1e-8)
        @test isapprox(corr.prefac, -4597019//80621568, atol = 1e-13)
    end

end