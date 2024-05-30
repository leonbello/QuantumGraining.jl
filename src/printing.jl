"""
printing.jl includes methods related to printing the different diagrams and terms.
"""

function Base.show(io::IO, d::DiagramNode)  
    write(io, "$(d.root) -> $(to_array(d))")
end
function Base.show(io::IO, d::NullNode)
    write(io, "$(d.root) -> NullNode")
end


"""
    symbolic_hamiltonian(gs::Vector, ops::Vector, Ω::Vector, t, τ)

Constructs a symbolic Hamiltonian for a quantum system.

# Arguments
- `gs::Vector`: Vector of symbols representing the coefficients of the Hamiltonian terms.
- `ops::Vector`: Vector of operators corresponding to each Hamiltonian term.
- `Ω::Vector`: Vector of frequencies for each Hamiltonian term.
- `t`: Time parameter.
- `τ`: Symbolic time parameter.

# Returns
- `terms`: Vector of symbolic Hamiltonian terms.

"""
function symbolic_hamiltonian(gs::Vector, ops::Vector, Ω::Vector, t, τ)
    terms = []
    for (g, op, ω) in zip(gs, ops, Ω)
        if isequal(ω, 0)
            ft = 1
        else
            ft = Symbolics.Term(exp, [im*ω*t])
        end

        term = to_symbol(g, τ)*ft
        if !isequal(term, 0)
            push!(terms, (term)*op)
        end
    end
    return terms
end

"""
    symbolic_hamiltonian(gs::Dict, Ω::Dict, t, τ)

Constructs a symbolic Hamiltonian for a quantum system.

# Arguments
- `gs::Dict`: A dictionary mapping operators to their corresponding coefficients.
- `Ω::Dict`: A dictionary mapping operators to their corresponding frequencies.
- `t`: The time parameter.
- `τ`: The time step parameter.

# Returns
A symbolic Hamiltonian for the quantum system.

"""
function symbolic_hamiltonian(gs::Dict, Ω::Dict, t, τ)
    gs_list = collect(values(gs))
    ops_list = collect(keys(gs))
    Ω_list = collect(values(Ω))

    return symbolic_hamiltonian(gs_list, ops_list, Ω_list, t, τ)
end

"""
    to_symbol(coeff::ContractionCoefficient, τ)

Give a symbolic representation of a `ContractionCoefficient` object.

# Arguments
- `coeff::ContractionCoefficient`: The `ContractionCoefficient` object to compute the symbol for.
- `τ`: The value of τ.

# Returns
Symbolic representation of the `ContractionCoefficient` object.

"""
function to_symbol(coeff::ContractionCoefficient, τ) #where {T <: Number}
    sym = 0
    for i in 1:length(coeff.prefacs)
        sym += to_symbol(Correction(coeff.prefacs[i], coeff.exponents[i], coeff.polys[i]), τ)
    end
    return sym
end

"""
    to_symbol(c::Correction, τ)

Give a symbolic representation of a `Correction` object.

# Arguments
- `c::Correction`: The `Correction` object to compute the symbol for.
- `τ`: The value of τ.

# Returns
Symbolic representation of the `Correction` object.

"""
function to_symbol(c::Correction, τ) #where {T <: Number}
    sym = c.prefac*exp(-0.5*τ^2*c.exponent)
    sym *= sum([isequal(c.poly[n], 0) ? 0 : c.poly[n]*(τ^(n-1)) for n in 1:c.order])
    return sym
end
