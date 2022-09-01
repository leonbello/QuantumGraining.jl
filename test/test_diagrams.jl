#using QuantumGraining
#using QuantumCumulants
module Tst 
    using Test                                                 # To make sure scopes don't get mixed up                             
    include("../src/diagrams.jl")
    include("../src/printing.jl")

    #= NullNode =#
    null = NullNode((3, 2))

    #= DiagramNode =#
    @test isequal(DiagramNode((3, 2)), DiagramNode((3, 2), (0, 0)))
    DiagramNode((3,2), (2, 1))

    node = DiagramNode((3, 2))
    node.left
    node.right
    child = DiagramNode(node.right)

    @test to_array(node.left)  == [(1, 0), (2, 2)]
    @test to_array(node.right) == [(0, 1), (3, 1)]
    
    # node_decomp() -- check by hand
    decomps = node_decomp(node)
    t = [ (el[1][1] + el[2][1], el[1][2] + el[2][2]) for el in decomps]
    @test all([el == (3, 2) for el in t])

    # get_diagrams() -- check that all sum up to the root
    node = DiagramNode((5, 4))
    diagrams = get_diagrams(node)
    @show diagrams
end



#= for automatic testing later
@testset "diagrams" begin
 
    
end #testset
=#