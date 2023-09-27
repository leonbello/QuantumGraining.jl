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
    include("../src/bvector.jl")
    include("../src/diagrams.jl")
    include("../src/bubble.jl")
    include("../src/diagram.jl")
    include("../src/poles.jl")
    include("../src/printing.jl")
    
   
    begin  #singularity in the first bubble should be omitted
        μ1 = UVec([0, 1, -2], special=true);  
        ν1 = DVec([])
        μ2 = UVec([1, 3]);
        ν2 = DVec([]);
        
        @show μ1[1]
        @show μ1[2]
        @show μ1[3]
        @show ν1
        @show μ1

        @show vec_factorial(μ1)
        @show vec_factorial(μ2)
        
        @show μ1.special
        @show μ1.type
        @show ν1.type
        @show μ1.poles
        @show ν1.poles

        @show sum(μ1)
        @show sum(ν1)
        @show sum(μ2)
        @show sum(ν2)

        @show length(μ1)
        @show length(ν1)

        @show μ1[2]
        @show size(μ1)

        @show sum(μ1)
        @show sum(ν1)
    end

    begin
        μ1 = UVec([0, 1, -1]; special=true);    
        ν1 = DVec([1]);

        @show vec_factorial(μ1)
        @show vec_factorial(ν1)
    end

    begin
        @cnumbers ω1 ω2 ω3 ω4
        μ1 = UVec([ω1, ω2, ω3, ω4]; special=true);
        
        result = vec_factorial(μ1)
        @show result
    end
end