using QuantumGraining
using QuantumCumulants
using Test

# DiagramNode
DiagramNode((3, 2)) == DiagramNode((3, 2), (3, 2))
DiagramNode((3,1), (2, 1))

node = DiagramNode((3, 2))
# right and left are mixed up
node.left
node.right


## Broken from here on...

# node_decomp()
a = node_decomp(node)
#=
@testset "diagrams" begin
 
    
end #testset
=#