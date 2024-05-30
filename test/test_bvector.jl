using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
using Test
using QuantumGraining

@testset "corrections" begin   
    begin  #singularity in the first bubble should be omitted
        μ1 = UVec([0, 1, -2], special=true);  
        ν1 = DVec([])
        μ2 = UVec([1, 3]);
        ν2 = DVec([]);
        
    end

    begin
        μ1 = UVec([0, 1, -1]; special=true);    
        ν1 = DVec([1]);
    end

    begin
        @cnumbers ω1 ω2 ω3 ω4
        μ1 = UVec([ω1, ω2, ω3, ω4]; special=true);
    end
end