module QuantumGraining

import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
#import TermInterface
#import SciMLBase

import ModelingToolkit
const MTK = ModelingToolkit

export DiagramNode, LeafNode,
        @definemodes, calculate_coeff,
        node_decomp, get_diagrams,
        coeff                                   # helper functions -- remember to remove after testing

export Contraction, 
        @definemodes, calculate_coeff,
        _Ω, coeff, _maxmodes                    # helper functions -- remember to remove after testing

include("diagrams.jl")
include("contractions.jl")

end
