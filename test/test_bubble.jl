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
    include("../src/expressions.jl")
    
   ## Bubble struct tests ##
   # Frequency constructor
    let
        shape = (3, 2)
        freqs = [1, 2, 3, 4, 5]
        bubble = Bubble(freqs, shape)
        @show bubble[1]
        @show bubble[2]
        @show bubble[3]
        @show bubble[4]
        @show bubble[5]
        @show bubble
        @show bubble.shape
        @show bubble.freqs
        @show bubble.up
        @show bubble.down
    end

    # BVector constructor
    let 
        μ1 = UVec([1, 2])
        ν1 = DVec([3, 4, 5])
        bubble = Bubble(μ1, ν1)
        @show bubble[1]
        @show bubble[2]
        @show bubble[3]
        @show bubble[4]
        @show bubble[5]
        @show bubble
        @show bubble.shape
        @show bubble.freqs
        @show bubble.up
        @show bubble.down
    end

    let 
        μ1 = [1, 2]
        ν1 = [3, 4, 5]
        bubble = Bubble(μ1, ν1)
        @show bubble[1]
        @show bubble[2]
        @show bubble[3]
        @show bubble[4]
        @show bubble[5]
        @show bubble
        @show bubble.shape
        @show bubble.freqs
        @show bubble.up
        @show bubble.down

        push!(μ1, 6)
        @show bubble
    end

    μ1 = [1, 2]
    ν1 = [3, 4, 5]
    bubble = Bubble(μ1, ν1; special=true)
    @show bubble
    @show bubble.special
    @show bubble.up.special
    @show bubble.down.special
end