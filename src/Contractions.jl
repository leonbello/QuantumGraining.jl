#Contractions.jl contains the diagramNode struct, 
#and getNodeDecompositions() function that gives all possible 
#decompositions of one node into two. 

struct diagramNode
    root::Tuple{Int64, Int64}
    val::Tuple{Int64, Int64}
    right::Any
    left::Any
    diagramNode(root, val) = new(root, val,0,0)

    function diagramNode(root, val)
        (val[1] > 0) ? left = diagramNode(root, (Int(val[1]-1), Int(val[2]))) : left = 0
        (val[2] > 0) ? right = diagramNode(root, (Int(val[1]), Int(val[2]-1))) : right = 0
        new(root, val, left, right)
        
    end
    
end


function getNodeDecompositions(breakNode::diagramNode)
    nodeDecomp_list = Array{}[]
    function NodeDecompositions(node)
        if (typeof(node) == diagramNode)
            nodeDecomp = [node.val, (node.root[1]-node.val[1], node.root[2]-node.val[2])]
            ((nodeDecomp in nodeDecomp_list) == false)  && (nodeDecomp[2][1] > 0) &&  push!(nodeDecomp_list, nodeDecomp)
            NodeDecompositions(node.left)
            NodeDecompositions(node.right)
        else
            return nodeDecomp_list
        #COMPLETE 
        end
    end
    return NodeDecompositions(breakNode)
end

