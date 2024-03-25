using Revise
using IterTools
using Symbolics
using Test
using QuantumGraining

@testset "contraction_coefficients" begin
    begin
        @variables μ τ
        c10 = contraction_coeff(1, 0, [μ])
        @show c10
        @test isequal(c10.exponents, [μ^2])
        @test isequal(c10.prefacs, [1])
    end

    begin
        @variables μ1 μ2 τ
        c20 = contraction_coeff(2, 0, [μ1, μ2])
        c11 = contraction_coeff(1, 1, [μ1, μ2])
        @show c20
        @show c11
        @test isequal(c20.exponents, [(μ1 + μ2)^2, (μ1^2 + μ2^2)])
        @test isequal(c20.prefacs, [-1/μ2, 1/μ2])
        @test isequal(c11.exponents, expand.(c20.exponents))
        @test isequal(simplify.(-c11.prefacs), c20.prefacs)
    end
end