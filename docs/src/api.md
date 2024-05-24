# [API](@id API)


## [BVector](@id API: BVector)

```@docs
DVec
UVec
BVector
```

## [Diagrams](@id API: Diagrams)

```@docs
DiagramNode
node_decomp
get_diagrams
Diagram
```

## [Poles](@id API: Poles)

```@docs
Contraction
count_poles
find_poles
find_all_poles
find_integer_solutions
reshape_sols
```


## [Contractions](@id API: Contractions)


```@docs
diagram_correction
contraction_coeff
calc_pole_corrections
Correction
ContractionCoefficient
split_freqs_into_bubbles
to_qc_symbol
```

## [Corrections](@id API: Corrections)
```@docs
merge_duplicate_exponents
simplify_contraction
```

## [Bubble](@id API: Bubble)
```@docs
Bubble
calculate_bubble_factor
```

## [Lindblad](@id API: Lindblad)
```@docs
effective_hamiltonian_term
effective_dissipator_term
repeated_combinations
gaussian_to_cutoff
drop_high_freqs
effective_hamiltonian
effective_dissipator
```

## [Printing](@id API: Printing)
```@docs
symbolic_hamiltonian
to_symbol
```

## [Convert](@id API: Convert)
```@docs
convert_expressions
hamiltonian_function
normal_ordered_dictionary
qc_convert
qnumber_to_qop
contraction_to_function
lindblad_function
```

## [Ordering](@id API: Ordering)
```@docs
expand_operators
expand_operator
group_operators
```
