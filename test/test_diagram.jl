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
    include("../src/bubble.jl")
    include("../src/diagrams.jl")
    include("../src/diagram.jl")
    include("../src/poles.jl")
    include("../src/printing.jl")
    
    μ1 = [0, 1, -2]
    ν1 = []

    μ2 = [1, 3]
    ν2 = []

    b1 = Bubble(μ1, ν1; special=true)
    b2 = Bubble(μ2, ν2)
    bubbles = [b1, b2]
    d = Diagram(bubbles)

    
    @show d
    @show d.bubbles
    @show d.bubbles[1].up
    @show d.bubbles[2].up
    @show d.shape
    @show d.freqs
    @show d.up_poles
    @show d.down_poles
end