"""
printing.jl includes methods related to printing the different diagrams and terms.
"""
function Base.show(io::IO, d::DiagramNode)  
    write(io, "$(to_array(d))")
end