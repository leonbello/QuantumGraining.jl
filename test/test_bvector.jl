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
    include("../src/printing.jl")
    include("../src/expressions.jl")
    include("../src/bvector.jl")
    include("../src/contractions.jl")
    
    # write tests for all basic possible diagram structures
    ## no singularities ##

    # only up-bubbles
    begin   # [(3, 0), (2, 0)], singularity in the first bubble should be omitted
        μ1 = UVec([0, 1, -2], special=true);  
        ν1 = DVec([]);
        μ2 = UVec([1, 3]);
        ν2 = DVec([]);
        ω = Diagram([(μ1, ν1), (μ2, ν2)]);
        
        @show typeof(ω)
        @show ω.shape
        @show ω.up_poles
        @show ω.down_poles
        @show ω.num_poles
        diagram_correction(ω)     
    end
end