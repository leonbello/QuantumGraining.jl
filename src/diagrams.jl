#=
diagrams.jl contains all recursive functions that help produce all child diagrams from a given contraction.
=#

abstract type AbstractDiagramNode end

"""
    LeafNode <: AbstractDiagramNode

A terminating node in the tree, its value is always chosen such that the rightmost bubble is (1, 0).

# Fields
- `root`: Root value of the tree.
- `val`: Value of the tree -- constant and always (1, 0).
"""
struct LeafNode <: AbstractDiagramNode
    root::Tuple{Int64, Int64}
    val::Tuple{Int64, Int64}
    
    LeafNode(root) = new(root, (root[1] - 1, root[2]))  # Value is chosen such that the rightmost bubble is (1,0)
end

"""
    DiagramNode

The DiagramNode struct is a type for one bubble in the diagram. 
It is defined as a recursive tree where the elements represent the different decompositions of one node into two. 

# Fields:
- `root::Tuple{Int, Int}`: The root of the tree, the largest common bubble defining the contraction.
- `val::Tuple{Int, Int}`: The value of the current node, determines how much of the bubble was broken to the left.
- `right::DiagramNode`: Pointer to the next right node of the tree, a diagram with one mode broken to the down-bubble.
- `left::DiagramNode`: Pointer to the next left node of the tree, a diagram with one mode broken to the up-bubble.
"""
struct DiagramNode <: AbstractDiagramNode
    root::Tuple{Int64, Int64}
    val::Tuple{Int64, Int64}
    left::AbstractDiagramNode
    right::AbstractDiagramNode

    diagramNode(root, val) = new(root, val, LeafNode(root), LeafNode(root))

    function DiagramNode(root, val)
        if ((root[1] < val[1]) || (root[2] < val[2]))
            throw(DomainError(root, "Root node must be smaller than node value!"))
        else
            left = (val[1] > root[1] - 1) ? DiagramNode(root, (val[1] - 1, val[2])) : LeafNode(root) # condition here is wrong
            right = (val[2] > root[2]) ? DiagramNode(root, (val[1], val[2] - 1)) : LeafNode(root)
            new(root, val, left, right)
        end 
    end
end
function DiagramNode(root::Tuple{Int64, Int64})
    return DiagramNode(root, root)
end

"""
    to_array(node::DiagramNode)

Given a node, returns the diagram in array form where each entry corresponds to a different bubble.
"""
function to_array(node::AbstractDiagramNode)
    return [node.val, (node.root[1] - node.val[1], node.root[2] - node.val[2])]
end


"""
    node_decomp(break_node::DiagramNode)

Uses the DiagramNode structure to give one level of decompositions explicitly using a recursive function.
In other words, gives all ways one can break a bubble into two.

# Arguments
- `node::DiagramNode`: Node to break down.

# Returns
- `decomp_list::Array{Int}`: a list of all nodes in the tree. 

"""
function node_decomp(node::DiagramNode)
    decomp_list = Array{}[]
    function _node_decomp(node::DiagramNode)
        decomp = to_array(node)
        if (decomp ∉ decomp_list && decomp isa DiagramNode)                           # changed for readability
            push!(decomp_list, decomp)
            _node_decomp(node.left)
            _node_decomp(node.right)
        else
            return decomp_list
        end
    end
    return _node_decomp(node)
end 

"""
    get_diagrams(init_diagram::Array)

Uses the above `node_decomp()` to get all possible diagrams recursively.

# Arguments
- `init_diagram::Array`: The initial diagram including the first level of breakdowns (two bubble diagrams)

# Returns
- A list of all possible diagrams for a given contractions.
"""
function get_diagrams(init_diagram::Array)
    #diagrams_list = Array{}[]
    diagrams_list = [init_diagram]
    function _get_diagrams(diagram::Array)
        if (last(diagram) == (1,0))                                         # stop rule: the stem-mode cannot be broken further.
            return diagrams_list
        else
            right_bubble = DiagramNode(last(diagram), last(diagram))        # create a new tree with the rightmost bubble as the root
            right_decomps = node_decomp(right_bubble)                       # returns all 2-bubble decompositions
            level_list = Array{}[]                                          # list of diagrams to be broken further
            for right_decomp in right_decomps
                child_diagram = copy(diagram)
                pop!(child_diagram)                                         # removes last bubble from the diagram
                push!(child_diagram, right_decomp...)
                push!(diagrams_list, child_diagram)
                push!(level_list, child_diagram)
            end
            for child_diagram in level_list
                _get_diagrams(child_diagram)
            end
        end
        return diagrams_list
    end
    #pushfirst!(diagrams_list, init_diagram)
    return _get_diagrams(init_diagram)
end
get_diagrams(root::DiagramNode) = get_diagrams(to_array(root))

