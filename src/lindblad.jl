"""
Lindblad.jl contains the functionalities to generate the final Lindbladian in operator form.
"""

"""
    effective_hamiltonian(c::Tuple{int, int})  
    
        Given a contraction, returns all contributing terms.
"""
using Combinatorics

function effective_hamiltonian(k::Int, ω::Array, h::Array, fmt=:QuantumCumulants)
    #Functionality to sum over all contractions of upto kth order 
    #Can just get the contractions we care about since we know what they should look like 
    #Assume omega and h are correctly ordered
    #Generate the expression - we assume calculate_coeff can give us the numerical coefficient
    #for a given set of frequencies and the corresponding diagram 

    #Write assumptions about contractions.jl
    #w and h are arrays of symbolic variables
    #Use tuple for w and h

    freq_list = collect(combinations(ω,k))
    h_list = collect(combinations(h,k))
    
    effective_ham_diagrams = get_diagrams(DiagramNode((k,0)))


    effective_hamiltonian = 0
    for j in range(1, length(freq_list), length(freq_list))
        for diagram in effective_ham_diagrams
            effective_hamiltonian += ((calculate_coeff(diagram, freq_list[j]) + calculate_coeff(diagram, reverse(-freq_list[j])))/2)*prod(h_list[j])
        end
    end
    return effective_hamiltonian
end


"""
function checkHam(c::Tuple{Int64, Int64})
    #Maybe a helper function 
    contraction = Contraction(c)
    effective_ham_diagrams = []
    for diagram in contraction.diagrams
        diagramInHam = 0
        for node in diagram
            if (node[2] == 0)
                diagramInHam = 1
            end
        end
        if (diagramInHam == 0)
            push!(diagram,effective_ham_diagrams)
        end
    end
end
"""
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
