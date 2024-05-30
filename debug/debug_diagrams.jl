using QuantumGraining

module Tst 
    using Test                                                 # To make sure scopes don't get mixed up                             
    
    #= DiagramNode =#
    @test isequal(DiagramNode((3, 2)), DiagramNode((3, 2), (0, 0)))
    DiagramNode((3,2), (2, 1))

    node = DiagramNode((3, 2))
    node.left
    node.right
    child = DiagramNode(node.right)

    # node_decomp() -- check by hand
    decomps = node_decomp(node)
    t = [ (el[1][1] + el[2][1], el[1][2] + el[2][2]) for el in decomps]
    @test all([el == (3, 2) for el in t])

    # get_diagrams() -- check that all sum up to the root
    node = DiagramNode((5, 4))
    diagrams = get_diagrams(node)
    @show typeof(diagrams)
    @show diagrams[1]

    diagram = diagrams[1]
    for (i, bubble) in enumerate(diagram)
        print(i)
        print(bubble)
    end
end

# end testset