"""
printing.jl includes methods related to printing the different diagrams and terms.
"""
function Base.show(io::IO, d::DiagramNode)  
    write(io, "$(d.root) -> $(to_array(d))")
end
function Base.show(io::IO, d::NullNode)
    write(io, "$(d.root) -> NullNode")
end