using Revise
using IterTools
using Symbolics
using Test
using QuantumGraining

@testset "contraction_coefficients" begin

    begin
        @variables μ1 μ2 τ
        c1 = contraction_coeff(2, 0, [μ1, μ2])
        c2 = c1 + 2*c1
        c3 = 3*c1      
        @test isequal(c2.exponents, merge_duplicate_exponents(c2).exponents)
        @test isequal(c2.exponents, c1.exponents)
        @test isequal(c3.prefacs, 3*c1.prefacs)
        @test isequal(c2.prefacs, 3*c1.prefacs)
    end

    begin
        @variables μ1 μ2 τ
        c20 = contraction_coeff(2, 0, [μ1, μ2])
        c10 = contraction_coeff(1, 0, [μ1])

        c3 = c20 + 2*c10
        @test issetequal(c3.exponents, unique([c10.exponents..., c20.exponents...]))
        @test isequal(c20.exponents, contraction_coeff(2, 0, [μ1, μ2]).exponents)
        @test isequal(c20.prefacs, contraction_coeff(2, 0, [μ1, μ2]).prefacs)
        @test issetequal(c10.exponents, contraction_coeff(1, 0, [μ1]).exponents)
        @test issetequal(c10.prefacs, contraction_coeff(1, 0, [μ1]).prefacs)
    end

    begin
        @variables μ τ
        c10 = contraction_coeff(1, 0, [μ])
        @test isequal(c10.exponents, [μ^2])
        @test isequal(c10.prefacs, [1])
    end

    begin
        @variables μ1 μ2 τ
        c20 = contraction_coeff(2, 0, [μ1, μ2])
        c11 = contraction_coeff(1, 1, [μ1, μ2])
        @test isequal(c20.exponents, [(μ1 + μ2)^2, (μ1^2 + μ2^2)])
        @test isequal(c20.prefacs, [-1/μ2, 1/μ2])
        @test isequal(c11.exponents, expand.(c20.exponents))
        @test isequal(simplify.(-c11.prefacs), c20.prefacs)
    end

    begin
        @variables μ1 μ2 τ
        c1 = contraction_coeff(2, 0, [μ1, μ2])
        c2 = contraction_coeff(1, 0, [μ1])
        c3 = c1 + c2
        @test issetequal(c3.exponents, unique([c1.exponents..., c2.exponents...]))
    end
end