using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining


@testset "dissipators" begin
    @variables g ωc ωa t τ
    Ω = [ωc + ωa, - ωc - ωa, ωc - ωa, - ωc + ωa]
    gvec = (g//2).*[1, 1, 1, 1]
    
    # Hilbert space definitions
    h_cav = FockSpace(:cavity)
    h_atom = NLevelSpace(:atom, (:g,:e))
    h = tensor(h_cav, h_atom)
    
    # Operator definitions
    @qnumbers a::Destroy(h) σ::Transition(h)
    σm = σ(:g, :e)
    σp = σ(:e, :g)
    σz = σ(:e, :e) - σ(:g, :g)
    σee = σ(:e, :e)
    hvec = [a*σm, a'*σp, a*σp, a'*σm]

    γ_eff_2, ω_eff_2 = effective_dissipator(hvec, gvec, Ω, 2)

    # (a*σm, a*σm):  2*ωa + 2*ωc
    @test isequal(ω_eff_2[(a*σm, a*σm)], 2*ωa + 2*ωc)
    @test isequal(length(γ_eff_2[(a*σm, a*σm)].exponents), 2)
    @test issetequal(γ_eff_2[(a*σm, a*σm)].exponents, [4((ωa + ωc)^2), 2((ωa + ωc)^2)])
    @test isequal(simplify(γ_eff_2[(a*σm, a*σm)].prefacs[1] - 1im*(g^2//2)*1/(ωa + ωc)), 0)
    @test isequal(simplify(γ_eff_2[(a*σm, a*σm)].prefacs[2] + 1im*(g^2//2)*1/(ωa + ωc)), 0)

    # ((a*σp), (a*σm)): 2ωc
    @test isequal(ω_eff_2[((a*σp), (a*σm))], 2*ωc)
    @test issetequal(γ_eff_2[((a*σp), (a*σm))].exponents, [(ωa + ωc)^2 + (-ωa + ωc)^2 + 2(ωa + ωc)*(-ωa + ωc), (ωa + ωc)^2 + (-ωa + ωc)^2])
    @test isequal(simplify(γ_eff_2[((a*σp), (a*σm))].prefacs[1] .- -1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a*σp), (a*σm))].prefacs[2] .- 1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test length(γ_eff_2[((a*σp), (a*σm))].exponents) == 2
    @test issetequal(γ_eff_2[((a*σp), (a*σm))].polys, [[1],[1]])

    # ((a'*σm), (a*σm)): 2ωa
    @test isequal(ω_eff_2[((a'*σm), (a*σm))], 2*ωa)
    @test issetequal(γ_eff_2[((a'*σm), (a*σm))].exponents, [(ωa + ωc)^2 + (ωa - ωc)^2 + 2(ωa + ωc)*(ωa - ωc), (ωa + ωc)^2 + (ωa - ωc)^2])
    @test length(γ_eff_2[((a'*σm), (a*σm))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σm), (a*σm))].prefacs[1] .- 1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a'*σm), (a*σm))].prefacs[2] .- -1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σm), (a*σm))].polys, [[1],[1]])

    # ((a'*σp), (a'*σp)): -2ωa - 2ωc
    @test isequal(ω_eff_2[((a'*σp), (a'*σp))], -2ωa - 2ωc)
    @test issetequal(γ_eff_2[((a'*σp), (a'*σp))].exponents, [4((-ωa - ωc)^2), 2((-ωa - ωc)^2)])
    @test length(γ_eff_2[((a'*σp), (a'*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σp), (a'*σp))].prefacs[1] .- -1im*(g^2//2)*1/(ωa + ωc)), 0)
    @test isequal(simplify((γ_eff_2[((a'*σp), (a'*σp))].prefacs[2] .- 1im*(g^2//2)*1/(ωa + ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σp), (a'*σp))].polys, [[1],[1]])

    # ((a*σp), (a'*σp)): -2ωa
    @test isequal(ω_eff_2[((a*σp), (a'*σp))], -2ωa)
    @test issetequal(γ_eff_2[((a*σp), (a'*σp))].exponents, [(-ωa - ωc)^2 + (-ωa + ωc)^2 + 2(-ωa + ωc)*(-ωa - ωc), (-ωa - ωc)^2 + (-ωa + ωc)^2])
    @test length(γ_eff_2[((a*σp), (a'*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a*σp), (a'*σp))].prefacs[1] .- -1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a*σp), (a'*σp))].prefacs[2] .- 1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a*σp), (a'*σp))].polys, [[1],[1]])

    # ((a'*σm), (a'*σp)): -2ωc
    @test isequal(ω_eff_2[((a'*σm), (a'*σp))], -2ωc)
    @test issetequal(γ_eff_2[((a'*σm), (a'*σp))].exponents, [(ωa - ωc)^2 + (-ωa - ωc)^2 + 2(ωa - ωc)*(-ωa - ωc), (ωa - ωc)^2 + (-ωa - ωc)^2])
    @test length(γ_eff_2[((a'*σm), (a'*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σm), (a'*σp))].prefacs[1] .- 1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a'*σm), (a'*σp))].prefacs[2] .- -1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σm), (a'*σp))].polys, [[1],[1]])

    # ((a*σm), (a*σp)): 2ωc
    @test isequal(ω_eff_2[((a*σm), (a*σp))], 2ωc)
    @test issetequal(γ_eff_2[((a*σm), (a*σp))].exponents, [(ωa + ωc)^2 + (-ωa + ωc)^2 + 2(ωa + ωc)*(-ωa + ωc), (ωa + ωc)^2 + (-ωa + ωc)^2])
    @test length(γ_eff_2[((a*σm), (a*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a*σm), (a*σp))].prefacs[1] .- -1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a*σm), (a*σp))].prefacs[2] .- 1im*(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a*σm), (a*σp))].polys, [[1],[1]])

    # ((a'*σp), (a*σp)): -2ωa
    @test isequal(ω_eff_2[((a'*σp), (a*σp))], -2ωa)
    @test issetequal(γ_eff_2[((a'*σp), (a*σp))].exponents, [(-ωa - ωc)^2 + (-ωa + ωc)^2 + 2(-ωa + ωc)*(-ωa - ωc), (-ωa - ωc)^2 + (-ωa + ωc)^2])
    @test length(γ_eff_2[((a'*σp), (a*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σp), (a*σp))].prefacs[1] .- -1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a'*σp), (a*σp))].prefacs[2] .- 1im*(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σp), (a*σp))].polys, [[1],[1]])


    # ((a*σp), (a*σp)): −2ωa + 2ωc
    @test isequal(ω_eff_2[((a*σp), (a*σp))], -2ωa + 2ωc)
    @test issetequal(γ_eff_2[((a*σp), (a*σp))].exponents, [4((-ωa + ωc)^2), 2((-ωa + ωc)^2)])
    @test length(γ_eff_2[((a*σp), (a*σp))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a*σp), (a*σp))].prefacs[1] .- -1im*(g^2//2)*1/(ωa - ωc)), 0)
    @test isequal(simplify(γ_eff_2[((a*σp), (a*σp))].prefacs[2] .- 1im*((g^2//2)*1/(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a*σp), (a*σp))].polys, [[1],[1]])

    # ((a*σm), (a'*σm)): 2ωa
    @test isequal(ω_eff_2[((a*σm), (a'*σm))], 2ωa)
    @test issetequal(γ_eff_2[((a*σm), (a'*σm))].exponents, [(ωa + ωc)^2 + (ωa - ωc)^2 + 2(ωa + ωc)*(ωa - ωc), (ωa + ωc)^2 + (ωa - ωc)^2])
    @test length(γ_eff_2[((a*σm), (a'*σm))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a*σm), (a'*σm))].prefacs[1] .- 1im*(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a*σm), (a'*σm))].prefacs[2] .- -1im*(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))), 0)
    @test issetequal(γ_eff_2[((a*σm), (a'*σm))].polys, [[1],[1]])

    # ((a'*σp), (a'*σm)): -2ωc
    @test isequal(ω_eff_2[((a'*σp), (a'*σm))], -2ωc)
    @test issetequal(γ_eff_2[((a'*σp), (a'*σm))].exponents, [(ωa - ωc)^2 + (-ωa - ωc)^2 + 2(ωa - ωc)*(-ωa - ωc), (ωa - ωc)^2 + (-ωa - ωc)^2])
    @test length(γ_eff_2[((a'*σp), (a'*σm))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σp), (a'*σm))].prefacs[1] .- 1im*(g^2//2)*ωc/((ωa - ωc)*(ωa + ωc))), 0)
    @test isequal(simplify(γ_eff_2[((a'*σp), (a'*σm))].prefacs[2] .- -1im*(g^2//2)*ωc/((ωa - ωc)*(ωa + ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σp), (a'*σm))].polys, [[1],[1]])

    # ((a'*σm), (a'*σm)): 2ωa−2ωc
    @test isequal(ω_eff_2[((a'*σm), (a'*σm))], 2ωa-2ωc)
    @test issetequal(γ_eff_2[((a'*σm), (a'*σm))].exponents, [4((ωa - ωc)^2), 2((ωa - ωc)^2)])
    @test length(γ_eff_2[((a'*σm), (a'*σm))].exponents) == 2
    @test isequal(simplify(γ_eff_2[((a'*σm), (a'*σm))].prefacs[1] .- 1im*(g^2//2)*1/(ωa - ωc)), 0)
    @test isequal(simplify(γ_eff_2[((a'*σm), (a'*σm))].prefacs[2] .- -1im*((g^2//2)*1/(ωa - ωc))), 0)
    @test issetequal(γ_eff_2[((a'*σm), (a'*σm))].polys, [[1],[1]])
end