using Revise
using QuantumCumulants
using SymbolicUtils
using IterTools
using Symbolics
using Test
using QuantumGraining

#@testset "contractions" begin
module Tst
    using Test
    using IterTools
    using QuantumGraining


    e1 = [1]
    pf1 = [1]
    c1 = ContractionCoefficient(e1, pf1)
    cr = -Correction(1, 1)


    e2 = [2, 3]
    pf2 = [2, 1]
    c2 = ContractionCoefficient(e2, pf2)

    c3 = c1 + c2
    c4 = 2*c3

    c5 = c4 - cr
end