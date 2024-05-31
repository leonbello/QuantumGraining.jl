using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining

@testset "lindblad" begin
    # Operator and space definitions
    begin
        @variables g ωc ωa
        #Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
        Ω = [ωc + ωa, - ωc - ωa, ωc - ωa, - ωc + ωa]
        gvec = (g/2).*[1, 1, 1, 1]

        h_cav = FockSpace(:cavity)
        h_atom = NLevelSpace(:atom, (:g,:e))
        h = tensor(h_cav, h_atom)

        @qnumbers a::Destroy(h) σ::Transition(h)
        σm = σ(:g, :e)
        σp = σ(:e, :g)
        σz = σ(:e, :e) - σ(:g, :g)
        σee = σ(:e, :e)
        hvec = [a*σm, a'*σp, a*σp, a'*σm]

        Σ = ωa + ωc
        Δ = ωc - ωa
    end

    # first-order Rabi-model - no RWA
    g_eff_1, Ω_eff_1 = effective_hamiltonian(hvec, gvec, Ω, 1; as_dict=true)

    @test issetequal(values(Ω_eff_1), [ωc + ωa, ωc - ωa, ωa - ωc, -ωc - ωa])
    @test issetequal(keys(g_eff_1), hvec)


    # a*σm: -ωa - ωc
    begin
        @test isequal(g_eff_1[a*σm].prefacs[1], g//2)
        @test isequal(expand(g_eff_1[a*σm].exponents[1]), expand((-ωa - ωc)^2))
        @test isequal(Ω_eff_1[a*σm], ωa + ωc)
    end

    # a'*σp: ωc + ωa
    begin
        @test isequal(g_eff_1[a'*σp].prefacs[1], g//2)
        @test isequal(expand(g_eff_1[a'*σp].exponents[1]), expand((ωc + ωa)^2))
        @test isequal(Ω_eff_1[a'*σp], -ωc - ωa)
    end

    # a*σp: ωa - ωc
    begin
        @test isequal(g_eff_1[a*σp].prefacs[1], g//2)
        @test isequal(expand(g_eff_1[a*σp].exponents[1]), expand((ωa - ωc)^2))
        @test isequal(Ω_eff_1[a*σp], -ωa + ωc)
    end

    # a'*σm: ωc - ωa
    begin
        @test isequal(g_eff_1[a'*σm].prefacs[1], g//2)
        @test isequal(expand(g_eff_1[a'*σm].exponents[1]), expand((ωc - ωa)^2))
        @test isequal(Ω_eff_1[a'*σm], -ωc + ωa)
    end

    ### RWA by dropping terms
    begin
        fsubs = Dict(
            ωa => 1,
            ωc => 1.01
        )
        g_loew_1, Ω_low_1 = drop_high_freqs(g_eff_1, Ω_eff_1, fsubs)
        @test issetequal(values(Ω_low_1), [ωc - ωa, ωa - ωc])
        @show collect(keys(Ω_low_1))

        g_low_1 = gaussian_to_cutoff(g_eff_1, Ω_low_1, fsubs; keep_small_exponents=true)
    end

    # a*σp: ωa - ωc
    begin
        @test isequal(g_low_1[a*σp].prefacs[1], g//2)
        @test isequal(expand(g_low_1[a*σp].exponents[1]), expand((ωa - ωc)^2))
    end

    # a'*σm: ωc - ωa
    begin
        @test isequal(g_low_1[a'*σm].prefacs[1], g//2)
        @test isequal(expand(g_low_1[a'*σm].exponents[1]), expand((ωc - ωa)^2))
    end
    
    # Second-order Rabi-model -- Hamiltonian terms
    g_eff_2, Ω_eff_2 = effective_hamiltonian(hvec, gvec, Ω, 2; as_dict=true)    
    @show keys(Ω_eff_2)    

    g1 = g^2/4*contraction_coeff(2, 0, [-ωa - ωc, ωa + ωc])
    g2 = g^2/4*contraction_coeff(2, 0, [ωa - ωc, ωc - ωa])
    begin
        # @show g1.prefacs
        # @show g^2/4*1/Σ .* [1, -1]
        # @show simplify.(g1.prefacs .- g^2/4*1/Σ .* [1, -1])
        @test issetequal(simplify.(g1.prefacs .- g^2/4*1/Σ.* [1, -1]), [0, 0])
    end

    begin
        @test issetequal(simplify.(g2.prefacs .- g^2/4*1/Δ.* [1, -1]), [0, 0])
    end

    # a'*a: 0
    begin
        gada = -(g1 - g2)
        @test isequal(simplify(g_eff_2[a'*a].exponents[1] - gada.exponents[1]), 0)
        @test isequal(simplify(g_eff_2[a'*a].exponents[2] - gada.exponents[2]), 0)
        @test isequal(simplify(g_eff_2[a'*a].exponents[3] - gada.exponents[3]), 0)

        @test isequal(simplify(g_eff_2[a'*a].prefacs[1] - gada.prefacs[1]), 0)
        @test isequal(simplify(g_eff_2[a'*a].prefacs[2] - gada.prefacs[2]), 0)
        @test isequal(simplify(g_eff_2[a'*a].prefacs[3] - gada.prefacs[3]), 0)
    end

    # σee: 0
    begin
        gσee = (g1 - g2)
        @test isequal(simplify(g_eff_2[σee].exponents[1] - gσee.exponents[1]), 0)
        @test isequal(simplify(g_eff_2[σee].exponents[2] - gσee.exponents[2]), 0)
        @test isequal(simplify(g_eff_2[σee].exponents[3] - gσee.exponents[3]), 0)

        @test isequal(simplify(g_eff_2[σee].prefacs[1] - gσee.prefacs[1]), 0)
        @test isequal(simplify(g_eff_2[σee].prefacs[2] - gσee.prefacs[2]), 0)
        @test isequal(simplify(g_eff_2[σee].prefacs[3] - gσee.prefacs[3]), 0)
    end

    # a'*a*σee: 0
    begin
        gadaσee = 2*(g1 - g2)
        @test isequal(simplify(g_eff_2[a'*a*σee].exponents[1] - gadaσee.exponents[1]), 0)
        @test isequal(simplify(g_eff_2[a'*a*σee].exponents[2] - gadaσee.exponents[2]), 0)
        @test isequal(simplify(g_eff_2[a'*a*σee].exponents[3] - gadaσee.exponents[3]), 0)

        @test isequal(simplify(g_eff_2[a'*a*σee].prefacs[1] - gadaσee.prefacs[1]), 0)
        @test isequal(simplify(g_eff_2[a'*a*σee].prefacs[2] - gadaσee.prefacs[2]), 0)
        @test isequal(simplify(g_eff_2[a'*a*σee].prefacs[3] - gadaσee.prefacs[3]), 0)
    end

    ### Third-order corrections
    g_eff_3, Ω_eff_3 = effective_hamiltonian(hvec, gvec, Ω, 3; as_dict=true)

    # (a'*a'*a'*σp): −ωa−3ωc
    test_exponents = [(-ωa - 3ωc)^2, (-ωa - ωc)^2 + (4//1)*(ωc^2), (ωa - ωc)^2 + 2((-ωa - ωc)^2)]
    test_prefacs = [0.271847, -1.33092, 1.05907]

    # The following results are tested with ωa=7/17, ωc=11/13, and g = 5/3
    begin
        @test isequal(Ω_eff_3[(a'*a'*a'*σp)], -ωa - 3*ωc)
    
        @test isequal(length(g_eff_3[(a'*a'*a'*σp)].exponents), 3)
    
        prefacs_num = [float.(substitute(g_eff_3[(a'*a'*a'*σp)].prefacs[i], Dict(ωa => 7/17, ωc => 11/13, g => 5/3))) for i in 1:3]
        @show g_eff_3[(a'*a'*a'*σp)].polys
        for i in 1:3
            @test isequal(g_eff_3[(a'*a'*a'*σp)].exponents[i], test_exponents[i])
            @test isapprox(prefacs_num[i], test_prefacs[i], rtol=0.001)
        end

        # the first element of polys should be 1
    end


    ### Fourth-order corrections
    g_eff_4, Ω_eff_4 = effective_hamiltonian(hvec, gvec, Ω, 4; as_dict=true)

    # (σee): 0
    test_exponents = [0, (ωa + ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + ωc)^2, 2*((ωa + ωc)^2 + (-ωa - ωc)^2), (ωa - ωc)^2 + (-ωa - ωc)^2 + (4//1)*(ωc^2), (ωa + ωc)^2 + (ωa - ωc)^2 + (-ωa - ωc)^2 + (-ωa + ωc)^2, 2*((ωa - ωc)^2 + (-ωa + ωc)^2)]
    test_prefacs = [[1.934400, 0, 0], [-0.912239, 0, -1.453626], [4.619092, 0, -1.453626], [0.242280], [-2.660242], [2.660242], [-5.883532]]
    # The following results are tested with ωa=7/17, ωc=11/13, and g = 5/3
    begin
        @test isequal(Ω_eff_4[σee], 0)

        @test isequal(length(g_eff_4[σee].exponents), 7)

        prefacs_num = [float.(substitute(g_eff_4[σee].prefacs[i]*g_eff_4[σee].polys[i], Dict(ωa => 0.41176470588235, ωc => 0.84615384615385, g => 1.66666666666667, g^2 => 2.77777777777778, g^4 => 7.71604938271605))) for i in 1:7]
        for i in 1:7
            @test isequal(g_eff_4[σee].exponents[i], test_exponents[i])
            for j in 1:length(test_prefacs[i])
                @test isapprox(prefacs_num[i][j], test_prefacs[i][j], rtol=0.001)
            end
        end
    end

    # (a*a*σee): 2ωc
    test_exponents = [(4//1)*(ωc^2), (ωa + ωc)^2 + (-ωa + ωc)^2, (ωa + ωc)^2 + (-ωa - ωc)^2 + (4//1)*(ωc^2), 2((ωa + ωc)^2) + (-ωa - ωc)^2 + (-ωa + ωc)^2, (ωa + ωc)^2 + (ωa - ωc)^2 + 2((-ωa + ωc)^2), (ωa - ωc)^2 + (-ωa + ωc)^2 + (4//1)*(ωc^2), (ωa + 3ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + 3ωc)^2]
    test_prefacs = [[5.072564, 0, 5.714011], [-7.415268, 0, 1.353134], [2.793791], [-2.793791], [0.317457], [-0.317457],[-0.385316],[2.728020]]
    # The following results are tested with ωa=7/17, ωc=11/13, and g = 5/3
    begin
        @test isequal(Ω_eff_4[a*a*σee], 2ωc)

        @test isequal(length(g_eff_4[a*a*σee].exponents), 8)

        prefacs_num = [float.(substitute(g_eff_4[a*a*σee].prefacs[i]*g_eff_4[a*a*σee].polys[i], Dict(ωa => 0.41176470588235, ωc => 0.84615384615385, g => 1.66666666666667, g^2 => 2.77777777777778, g^4 => 7.71604938271605))) for i in 1:8]
        for i in 1:8
            @test isequal(g_eff_4[a*a*σee].exponents[i], test_exponents[i])
            for j in 1:length(test_prefacs[i])
                @test isapprox(prefacs_num[i][j], test_prefacs[i][j], rtol=0.001)
            end
        end
    end                     

    # (a'*a*σee): 0
    test_exponents = [0, (ωa + ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + ωc)^2, 2((ωa + ωc)^2 + (-ωa - ωc)^2), (ωa + ωc)^2 + (-ωa + ωc)^2 + (4//1)*(ωc^2), (ωa + ωc)^2 + (ωa - ωc)^2 + (-ωa - ωc)^2 + (-ωa + ωc)^2, 2((ωa - ωc)^2 + (-ωa + ωc)^2)]
    test_prefacs = [[9.83081956], [-2.544840, 0, -5.814503], [15.279030, 0, -5.814503], [0.9691197], [-10.640968], [10.640968],[-23.53412897]]
    # The following results are tested with ωa=7/17, ωc=11/13, and g = 5/3
    begin
        @test isequal(Ω_eff_4[a'*a*σee], 0)

        @test isequal(length(g_eff_4[a'*a*σee].exponents), 7)

        prefacs_num = [float.(substitute(g_eff_4[a'*a*σee].prefacs[i]*g_eff_4[a'*a*σee].polys[i], Dict(ωa => 0.41176470588235, ωc => 0.84615384615385, g => 1.66666666666667, g^2 => 2.77777777777778, g^4 => 7.71604938271605))) for i in 1:7]
        for i in 1:7
            @test isequal(g_eff_4[a'*a*σee].exponents[i], test_exponents[i])
            for j in 1:length(test_prefacs[i])
                @test isapprox(prefacs_num[i][j], test_prefacs[i][j], rtol=0.001)
            end
        end
    end

    # (a'*a'*a'*a'*σee): -4ωc
    test_exponents = [(16//1)*(ωc^2), (-ωa - ωc)^2 + (ωa - 3ωc)^2, (8//1)*(ωc^2), (ωa - ωc)^2 + (-ωa - 3ωc)^2, (ωa - ωc)^2 + (-ωa - ωc)^2 + (4//1)*(ωc^2), 2((ωa - ωc)^2 + (-ωa - ωc)^2)]
    test_prefacs = [-0.231682, -0.213043, 1.330121, 0.444725, -2.660242, 1.330121]
    
    # The following results are tested with ωa=7/17, ωc=11/13, and g = 5/3
    begin
        @test isequal(Ω_eff_4[(a'*a'*a'*a'*σee)], -4*ωc)
    
        @test isequal(length(g_eff_4[(a'*a'*a'*a'*σee)].exponents), 6)
    
        prefacs_num = [float.(substitute(g_eff_4[(a'*a'*a'*a'*σee)].prefacs[i], Dict(ωa => 7/17, ωc => 11/13, g => 5/3))) for i in 1:6]
        @show g_eff_4[(a'*a'*a'*a'*σee)].polys
        for i in 1:6
            @test isequal(g_eff_4[(a'*a'*a'*a'*σee)].exponents[i], test_exponents[i])
            @test isapprox(prefacs_num[i], test_prefacs[i], rtol=0.001)
        end

        # the first element of polys should be 1
    end
end