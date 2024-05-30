using Revise
using IterTools
using Symbolics
using Test
using QuantumGraining


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


left = 2
right = 3
freqs = ones(left + right)
#contraction_coeff(left, right, freqs::Array)

## Error in inputing the second frequency into the diagram -- for some reason poly is zero
ωs = [1, 0]
diagrams = get_diagrams(DiagramNode((2, 0)))
ω_test = split_freqs_into_bubbles(ωs, diagrams[2])
corr = diagram_correction(ω_test)
@show corr.poly