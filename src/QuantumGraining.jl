module QuantumGraining

# "we do things not because they are easy, but because we thought they were"

import SymbolicUtils
import SymbolicUtils: substitute
import QuantumCumulants
import Symbolics
import IterTools
#import SymPy
#import SciMLBase

import ModelingToolkit
const MTK = ModelingToolkit

# bvector.jl
export DVec, UVec, BVector

# decomp.jl
export DiagramNode, node_decomp!,
        node_decomp, get_diagrams, Diagram

# poles.jl
export Contraction, count_poles,
        find_poles, find_all_poles, find_integer_solutions,
        reshape_sols

# contractions.jl
export diagram_correction, contraction_coeff, pole_corrections,
        Correction, ContractionCoefficient,
        split_freqs_into_tuples, to_qc_symbol, simple_factors, expansion_factors, poly_multiplication

# corrections.jl
export merge_duplicate_exponents, simplify

# bubble.jl
export Bubble, bubble_factor

# lindblad.jl
export effective_hamiltonian_term, effective_dissipator_term, repeated_combinations,
        gaussian_to_cutoff, drop_high_freqs, effective_hamiltonian, effective_dissipator

# symbols.jl
export symbolic_hamiltonian, to_symbol

# convert.jl
export convert_expressions, hamiltonian_function, normal_ordered_dictionary, qc_convert, qnumber_to_qop,
        contraction_to_function, lindblad_function

# ordering.jl
export expand_operators, expand_operator, group_operators

include("bvector.jl")
include("bubble.jl")
include("diagram.jl")
include("corrections.jl")
include("contractions.jl")
include("decomp.jl")
include("lindblad.jl")
include("poles.jl")
include("symbols.jl")
include("utils.jl")
include("convert.jl")
include("ordering.jl")
end
