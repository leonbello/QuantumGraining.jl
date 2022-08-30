module QuantumGraining

import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
<<<<<<< HEAD
#import TermInterface

#import SciMLBase
=======
import TermInterface

import SciMLBase
>>>>>>> contraction-coeffecients

import ModelingToolkit
const MTK = ModelingToolkit

<<<<<<< HEAD
export DiagramNode, LeafNode,
        @definemodes, calculate_coeff,
        node_decomp, get_diagrams,
        coeff                                   # helper functions -- remember to remove after testing

export Contraction, 
        @definemodes, calculate_coeff,
        _Ω, coeff, _maxmodes                    # helper functions -- remember to remove after testing


include("Diagrams.jl")
include("Contractions.jl")

=======
export DiagramNode

include("Contractions.jl")
include("Diagrams.jl")
>>>>>>> contraction-coeffecients

end
