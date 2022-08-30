module QuantumGraining

import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
<<<<<<< Updated upstream
<<<<<<< HEAD
#import TermInterface

#import SciMLBase
=======
import TermInterface

import SciMLBase
>>>>>>> contraction-coeffecients
=======
#import TermInterface

#import SciMLBase
>>>>>>> Stashed changes

import ModelingToolkit
const MTK = ModelingToolkit

<<<<<<< Updated upstream
<<<<<<< HEAD
=======
>>>>>>> Stashed changes
export DiagramNode, LeafNode,
        @definemodes, calculate_coeff,
        node_decomp, get_diagrams,
        coeff                                   # helper functions -- remember to remove after testing

export Contraction, 
        @definemodes, calculate_coeff,
        _Ω, coeff, _maxmodes                    # helper functions -- remember to remove after testing
<<<<<<< Updated upstream


include("Diagrams.jl")
include("Contractions.jl")
=======
>>>>>>> Stashed changes

=======
export DiagramNode

include("Diagrams.jl")
<<<<<<< Updated upstream
>>>>>>> contraction-coeffecients
=======
include("Contractions.jl")

>>>>>>> Stashed changes

end
