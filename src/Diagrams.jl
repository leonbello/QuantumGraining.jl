#=
This file should contain the recursive functions that help produce 
all child diagrams from a given Contractions

The DiagramNode struct is a type for one bubble in the diagram 
It is defined as a recursive tree where the elements represent the different decompositions
of one node into two. 

getNodeDecompositions uses the DiagramNode structure to give the decompositions explicitly 

getDiagrams uses the above function to get all possible diagrams recursively
=#


# TODO:
# - right and left should be DiagramNode types.

struct DiagramNode
    root::Tuple{Int64, Int64}
    val::Tuple{Int64, Int64}
    right::Any
    left::Any
    DiagramNode(root, val) = new(root, val,0,0)

    function DiagramNode(root, val)
        if ((root[1] < val[1]) || (root[2] < val[2]))
            throw(DomainError(root, "Root node must be smaller than Val node"))
        else
            (val[1] > 0) ? left = DiagramNode(root, (Int(val[1]-1), Int(val[2]))) : left = 0
            (val[2] > 0) ? right = DiagramNode(root, (Int(val[1]), Int(val[2]-1))) : right = 0
            new(root, val, left, right)
        end
        
    end
    
end


function getNodeDecompositions(breakNode::DiagramNode)
    nodeDecomp_list = Array{}[]
    function NodeDecompositions(node)
        if (typeof(node) == DiagramNode) 
            nodeDecomp = [node.val, (node.root[1]-node.val[1], node.root[2]-node.val[2])]
            ((nodeDecomp in nodeDecomp_list) == false) && (node.val != (0,0)) && (nodeDecomp[2][1] > 0)  && push!(nodeDecomp_list, nodeDecomp)
            NodeDecompositions(node.left)
            NodeDecompositions(node.right)
        else
            return nodeDecomp_list
        end
    end
    return NodeDecompositions(breakNode)
end


function getDiagrams(init_diagram::Array)
    diagramsList = Array{}[]
    function getChildDiagrams(diagram::Array)
        if (last(diagram) == (1,0))
            return diagramsList
        else
            rightBubble = DiagramNode(last(diagram), last(diagram))
            rightDecompList = getNodeDecompositions(rightBubble)
            levelDiagramList = Array{}[]
            for rightDecomp in rightDecompList
                childDiagram = copy(diagram)
                pop!(childDiagram)
                push!(childDiagram, rightDecomp[1])
                push!(childDiagram, rightDecomp[2])
                push!(diagramsList, childDiagram)
                push!(levelDiagramList, childDiagram)
            end
            for childDiagram in levelDiagramList
                getChildDiagrams(childDiagram)
            end


        end
        return diagramsList
    end
    pushfirst!(diagramsList, init_diagram)
    return getChildDiagrams(init_diagram)
end
