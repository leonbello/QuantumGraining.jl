using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
using Test
using QuantumGraining

@testset "bubble" begin
    
   
    # Frequency constructor
    begin
        shape = (3, 2)
        freqs = [1, 2, 3, 4, 5]
        bubble = Bubble(freqs, shape)
        @test bubble.shape == (3, 2)
        @test bubble.freqs == freqs
        @test bubble.up == [1, 2, 3]
        @test bubble.down == [4, 5]
    end

    # BVector constructor
    begin
        μ1 = UVec([1, 2])
        ν1 = DVec([3, 4, 5])
        bubble = Bubble(μ1, ν1)
        @test bubble.shape == (2, 3)
        @test bubble.freqs == [1, 2, 5, 4, 3]
        @test bubble.up == [1, 2]
        @test bubble.down == [3, 4, 5]
        @test bubble.special == false
    end

    begin
        μ1 = [1, 2]
        ν1 = [3, 4, 5]
        bubble = Bubble(μ1, ν1; special=true)
        @test bubble.shape == (2,3)
        @test bubble.freqs == [1, 2, 5, 4, 3]
        @test bubble.special == true
    end
end