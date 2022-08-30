"""
Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
"""

"""
    effective_hamiltonian(c::Tuple{int, int})  
    
        Given a contraction, returns all contributing terms.
"""
function effective_hamiltonian(c::Tuple{int, int}, fmt=:QuantumCumulants)
    print("Hello")
end

"""
    effective_hamiltonian(d::Diagram) - Given a diagram object, returns all contributing terms.
    effective_hamiltonian(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""

"""
* effective_dissipator(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_dissipator(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_dissipator(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""

"""
* effective_lindblad(c::Tuple{int, int}) - Given a contraction, returns all contributing terms.
* effective_lindblad(d::Diagram) - Given a diagram object, returns all contributing terms.
* effective_lindblad(d::Array{Tuple{Int, Int}}) - Given a diagram in an array format, returns all contributing terms.
"""