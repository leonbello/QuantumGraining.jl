#=
diagrams.jl contains all recursive functions that help produce all child diagrams from a given contraction.
=#
abstract type AbstractDiagramNode end

"""
    NullNode <: AbstractDiagramNode

A terminating node in the tree, has no value. Here so we can check if node can be decomposed using `isa DiagramNode`

# Fields
- `root`: Root value of the tree.

"""
struct NullNode <: AbstractDiagramNode
    root::Tuple{Int64, Int64}
end

"""
    DiagramNode

The DiagramNode struct is a type for one bubble in the diagram. 
It is defined as a recursive tree where the elements represent the different decompositions of one node into two. 

# Fields:
- `root::Tuple{Int, Int}`: The root of the tree, the largest common bubble defining the contraction.
- `val::Tuple{Int, Int}`: The value of the current node, determines how much of the bubble was broken to the left.
- `rightmost::Tuple{Int, Int}`: The value of the rightmost bubble.
- `right::DiagramNode`: Pointer to the next right node of the tree, a diagram with one mode broken to the down-bubble.
- `left::DiagramNode`: Pointer to the next left node of the tree, a diagram with one mode broken to the up-bubble.
"""
struct DiagramNode <: AbstractDiagramNode
    root::Tuple{Int64, Int64}
    val::Tuple{Int64, Int64}
    rightmost::Tuple{Int64, Int64}
    left::AbstractDiagramNode
    right::AbstractDiagramNode

    function DiagramNode(root, val)
        if ((root[1] < val[1]) || (root[2] < val[2]))
            throw(DomainError(root, "Root node must be smaller than node value!"))
        else
            rightmost = (root[1] - val[1], root[2] - val[2])
            # the nodes (2, N) and (N, 1) are the last ones we can break.
            left  = (rightmost[1] >= 2) ? DiagramNode(root, (val[1] + 1, val[2])) : NullNode(root)
            right = (rightmost[2] >= 1) ? DiagramNode(root, (val[1], val[2] + 1)) : NullNode(root)
            new(root, val, rightmost, left, right)
        end 
    end
end
function DiagramNode(root::Tuple{Int64, Int64})
    return DiagramNode(root, (0, 0))
end
function DiagramNode(node::DiagramNode)
    DiagramNode(node.rightmost)
end

"""
    to_array(node::DiagramNode)

Given a node, returns the diagram in array form where each entry corresponds to a different bubble.
"""
function to_array(node::DiagramNode)
    return [node.val, node.rightmost]
end
to_array(node::NullNode) = []

"""
    node_decomp(node::DiagramNode)
    node_decomp!(node::AbstractDiagramNode, decomp_list)

Uses the DiagramNode structure to give one level of decompositions explicitly using a recursive function.
In other words, gives all ways one can break a bubble into two.

# Arguments
- `node::DiagramNode`: Node to break down.

# Returns
- `decomp_list::Array{Int}`: a list of all nodes in the tree. 
"""
function node_decomp!(node::AbstractDiagramNode, decomp_list)
    if (node isa DiagramNode)
        decomp = to_array(node)
        if (decomp ∉ decomp_list)                                        # changed for readability
            push!(decomp_list, decomp)
        end
        node_decomp!(node.left, decomp_list)
        node_decomp!(node.right, decomp_list)
    end
end

function node_decomp(node::DiagramNode)
    decomp_list = []
    node_decomp!(node, decomp_list)
    return decomp_list
end 

"""
    get_diagrams(init_diagram::Array)

Uses the above `node_decomp()` to get all possible diagrams recursively.

# Arguments
- `node::AbstractDiagramNode`: The initial diagram including the first level of breakdowns (two bubble diagrams)

# Returns
- A list of all possible diagrams for a given contractions.
"""
function get_diagrams(node::AbstractDiagramNode)
    diagrams_list = []
    if (node isa DiagramNode)
        node_decomp!(node, diagrams_list)                              # add all 2-bubble decompoisition to the list
        child_set = []
        reference_node = node
        current_node = node
        while true
            if current_node.left isa DiagramNode
                push!(child_set, current_node.left)
                current_node = current_node.left
            elseif reference_node.right isa DiagramNode
                push!(child_set, reference_node.right)
                current_node = reference_node.right
                reference_node = reference_node.right
            else
                break
            end
        end

        for child in child_set                                         # for all child nodes
            if (child isa DiagramNode)
                decomp = DiagramNode(child)                            # build a new decomposition tree out of the rightmost bubble
                child_list = get_diagrams(decomp)                      # get all diagrams of the child node
                for c in child_list
                    pushfirst!(c, child.val)                           # for each decomposition, add the left bubble
                end
                push!(diagrams_list, child_list...)
            end
        end
    end
    diagrams_list = pushfirst!(filter!(x-> ((0,0) ∉ x), diagrams_list),[node.root]) #Removing redundant diagrams with (0,0) terms
    return unique!(diagrams_list)
end

function get_diagrams(l::Int, r::Int)
    return get_diagrams(DiagramNode((l, r)))
end