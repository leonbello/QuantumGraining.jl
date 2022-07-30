module Contractions

mutable struct diagramNode
    root::Tuple{Int, Int}
    val::Tuple{Int, Int}
    right::diagramNode
    left::diagramNode

    function diagramNode(root, val)
        (val[1] > 0) ? left = diagramNode(root, (Int(val[1]-1), Int(val[2]))) : left = 0
        (val[2] > 0) ? right = diagramNode(root, (Int(val[1]), Int(val[2]-1))) : right = 0
        new(root, val, left, right)
        
    end

end

end