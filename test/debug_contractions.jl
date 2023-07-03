using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
using Test

include("../src/diagrams.jl")
include("../src/bvector.jl")
include("../src/bubble.jl")
include("../src/diagram.jl")
include("../src/poles.jl")
include("../src/contractions.jl")
include("../src/printing.jl")


# tests for calculate_bubble_factor
μ = [1, 2]
ν = [3, 4, 5]
b = Bubble(μ, ν; special=true)

num_poles = length(b.up.poles) + length(b.down.poles)
sols = find_integer_solutions(3, num_poles)    
sols = reshape_sols(sols, num_poles, 1) 
calculate_bubble_factor(b, sols, b.up.poles, b.down.poles)


# tests for diagram correction functions
μ1 = [0, 1, -2];  
ν1 = [];
μ2 = [1, 3];
ν2 = [];
d = [(μ1, ν1), (μ2, ν2)];
diagram_correction(d) 
