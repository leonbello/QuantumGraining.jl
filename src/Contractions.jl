module Contractions

mutable struct diagramNode
    root::AbstractArray{Tuple{Int, Int}, N}
    val::Tuple{Int, Int}
    right::diagramNode
    left::diagramNode

    function getDiagram(node::diagramNode)
        partitions = [] #Get all partitions of node.val
        diagram = !append(partitions, node.root)
        
    end

end

#Write recursive algorithm that does the partitions at each level 



end