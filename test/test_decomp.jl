using QuantumGraining

@testset "decomp" begin
    @test get_diagrams(2, 0) == [[(2, 0)], [(1, 0), (1, 0)]]
    @test get_diagrams(1, 1) == [[(1, 1)], [(0, 1), (1, 0)]]
end