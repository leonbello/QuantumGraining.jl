using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining


@testset "lindblad" begin
    # Operator and space definitions
    begin
        @variables g ωc ωa
        Ω = [-ωc - ωa, ωc + ωa, -ωc + ωa, ωc - ωa]
        gvec = (g/2).*[1, 1, 1, 1]

        h_cav = FockSpace(:cavity)
        h_atom = NLevelSpace(:atom, (:g,:e))
        h = tensor(h_cav, h_atom)

        @qnumbers a::Destroy(h) σ::Transition(h)
        σm = σ(:e, :g)
        σp = σ(:g, :e)
        hvec = [a*σm, a'*σp, a*σp, a'*σm]

        Σ = ωa + ωc
        Δ = ωc - ωa
    end

    # first-order Rabi-model - no RWA
    ops_eff, g_eff, Ω_eff = effective_hamiltonian(hvec, gvec, Ω, 1)

    begin
        @test issetequal(Ω_eff, [ωc + ωa, ωc - ωa, ωa - ωc, -ωc - ωa])
        @test isequal(g_eff[1].prefacs[1], g)
        @test isequal(g_eff[1].exponents[1], (-ωc - ωa)^2)
        @test isequal(ops_eff[1], a*σm)       
    end

    begin
        @test isequal(g_eff[2].prefacs[1], g)
        @test isequal(g_eff[2].exponents[1], (ωc + ωa)^2)
        @test isequal(ops_eff[2], a'*σp)
    end

    begin
        @test isequal(g_eff[3].prefacs[1], g)
        @test isequal(g_eff[3].exponents[1], (ωa - ωc)^2)        
        @test isequal(ops_eff[3], a*σp)
    end

    begin
        @test isequal(g_eff[4].prefacs[1], g)
        @test isequal(g_eff[4].exponents[1], (ωc - ωa)^2)
        @test isequal(ops_eff[4], a'*σm)
    end

    ### RWA by dropping terms
    begin
        fsubs = Dict(
            ωa => 1,
            ωc => 1.01
        )
        rwa, Ω_rwa = drop_high_freqs(Ω_eff, fsubs)
        @test issetequal(Ω_rwa, [ωc - ωa, ωa - ωc])
        
        ops_rwa = ops_eff[rwa]
        g_rwa = g_eff[rwa]
        g_cut = gaussian_to_cutoff(g_rwa, fsubs)  
    end

    begin
        @test isequal(g_rwa[1].prefacs[1], g)
        @test isequal(g_rwa[1].exponents[1], (ωa - ωc)^2)
        @test isequal(g_cut[1].exponents[1], 0)
    end

    begin
        @test isequal(g_rwa[2].prefacs[1], g)
        @test isequal(g_rwa[2].exponents[1], (ωc - ωa)^2)
        @test isequal(g_cut[2].exponents[1], 0)
    end
    
    # Second-order Rabi-model -- Hamiltonian terms
    ops_eff, g_eff, Ω_eff = effective_hamiltonian_term(hvec, gvec, Ω, 2)    

    # -ω, -ω: 
    begin        
        @test isequal(Ω_eff[1], -2*ωa - 2*ωc)
        @test isequal(ops_eff[1], 0) # equals 0 
    end

    # ω, -ω:
    begin    
        @test isequal(Ω_eff[2], 0)
        @test isequal(simplify(ops_eff[2] - a'*σp*a*σm), 0)
        @test issetequal(expand.(g_eff[2].exponents), [0, 2*ωa^2 + 4*ωa*ωc + 2*ωc^2])
        @test issetequal(simplify.(g_eff[2].prefacs .- -[-g^2//2*1/(ωa + ωc), g^2//2*1/(ωa + ωc)]), [0,0])
    end

    # -Δ, -ω
    begin
        @test isequal(Ω_eff[3], -(ωc - ωa) - (ωc + ωa))
        @test isequal(simplify(ops_eff[3] - a*σp*a*σm), 0)
        @test issetequal(
            expand.(g_eff[3].exponents),
            expand.([Δ^2 + Σ^2, (Σ + Δ)^2])
        )
        @test isequal(
            simplify(g_eff[3].prefacs[1] - g^2//4*(1/Σ - 1/Δ)),
            0
        )
        @test isequal(
            simplify(g_eff[3].prefacs[2] - -g^2//4*(1/Σ - 1/Δ)),
            0
        )
    end

    # -ω, Δ:
    begin
        @test isequal(Ω_eff[4], -ωa - ωc + ωc - ωa)
        @test isequal(simplify(ops_eff[4] - a*σm*a'*σm), 0) # equals 0       
    end
    
    # -ω, ω:    
    begin
        @test isequal(Ω_eff[5], 0)
        @test isequal(simplify(ops_eff[5] - a*σm*a'*σp), 0)
        @test issetequal(expand.(g_eff[5].exponents), [0, 2*ωa^2 + 4*ωa*ωc + 2*ωc^2])
        @test issetequal(simplify.(g_eff[5].prefacs .- [-g^2//2*1/(ωa + ωc), g^2//2*1/(ωa + ωc)]), [0,0])
    end

    # ω, ω:
    begin
        @test isequal(Ω_eff[6], 2ωa + 2ωc)
        @test isequal(ops_eff[6], 0)  # equals 0  
    end

    # ω, -Δ: 
    begin   
        @test isequal(Ω_eff[7], (ωa + ωc) - (ωc - ωa))
        @test isequal(ops_eff[7], a'*σp*a*σp) # equals 0
    end
    
    # Δ, ω:
    begin
        @test isequal(Ω_eff[8], (ωa + ωc) + (ωc - ωa))
        @test isequal(ops_eff[8], a'*σm*a'*σp)
        @test issetequal(
            expand.(g_eff[8].exponents),
            expand.([Δ^2 + Σ^2, (Σ + Δ)^2])
        )
        @test isequal(
            simplify(g_eff[8].prefacs[1] - -g^2//4*(1/Σ - 1/Δ)),
            0
        )
        @test isequal(
            simplify(g_eff[8].prefacs[2] - g^2//4*(1/Σ - 1/Δ)),
            0
        )
    end
    
    # -ω, -Δ:
    begin
        @test isequal(Ω_eff[9], -(ωc + ωa) - (ωc - ωa))
        @test isequal(ops_eff[9], a*σm*a*σp)
        @test issetequal(
            expand.(g_eff[9].exponents),
            expand.([Δ^2 + Σ^2, (Σ + Δ)^2])
        )
        @test isequal(
            simplify(g_eff[9].prefacs[1] - -g^2//4*(1/Σ - 1/Δ)),
            0
        )
        @test isequal(
            simplify(g_eff[9].prefacs[2] - g^2//4*(1/Σ - 1/Δ)),
            0
        )
    end

    # -Δ, ω:
    begin    
        @test isequal(Ω_eff[10], -(ωc - ωa) + (ωc + ωa))
        @test isequal(ops_eff[10], a*σp*a'*σp) # equals 0
    end

    # -Δ, -Δ:
    begin
        @test isequal(simplify(Ω_eff[11] - -2*(ωc - ωa)), 0)
        @test isequal(ops_eff[11], a*σp*a*σp)  # equals 0      
    end

    # Δ, Δ:
    begin
        @test isequal(simplify(Ω_eff[16] - 2*(ωc - ωa)), 0)
        @test isequal(ops_eff[16], a'*σm*a'*σm) # equals 0        
    end

    # Δ, -Δ:
    begin
        @test isequal(Ω_eff[12], -(ωc - ωa) + (ωc - ωa))
        @test isequal(ops_eff[12], a'*σm*a*σp)
        @test issetequal(
            simplify.(expand.(g_eff[12].exponents)), 
            simplify.(expand.([2*(ωc - ωa)^2, 0]))
        )
        @test issetequal(
            simplify.(g_eff[12].prefacs .- g^2//2*1/(ωc - ωa).*[1, -1]),
            0
        )         
    end

    # -Δ, Δ:
    begin
        @test isequal(Ω_eff[15], -(ωc - ωa) + (ωc - ωa))
        @test isequal(simplify(ops_eff[15] - a*σp*a'*σm), 0)
        @test issetequal(
            simplify.(expand.(g_eff[15].exponents)), 
            simplify.(expand.([2*(ωc - ωa)^2, 0]))
        )
        @test issetequal(
            simplify.(g_eff[15].prefacs .- -g^2//2*1/(ωc - ωa).*[1, -1]),
            [0, 0]
        )        
    end

    # Δ, -ω:
    begin    
        @test isequal(Ω_eff[13], (ωc - ωa) - (ωa + ωc))
        @test isequal(ops_eff[13], a'*σm*a*σm) # equals 0
    end

    # ω, Δ:
    begin                
        @test isequal(Ω_eff[14], (ωc - ωa) + (ωa + ωc))
        @test isequal(
            simplify(ops_eff[14] - a'*σp*a'*σm),
            0
        )
        @test issetequal(
            expand.(g_eff[14].exponents),
            expand.([(ωc + ωa)^2 + (ωc - ωa)^2, ((ωc + ωa) + (ωc - ωa))^2])
        )
        @test issetequal(
            simplify.(g_eff[14].prefacs .- g^2//4*(1/(ωc - ωa) - 1/(ωc + ωa)).*[-1, 1]),
            [0, 0]
        )
    end

    # Second-order Rabi-model -- Dissipators
    begin
        diss_eff, γ_eff, Ω_eff = effective_dissipator_term(hvec, gvec, Ω, 2)
        rwa, Ω_rwa = drop_high_freqs(Ω_eff, fsubs)
        @show Ω_rwa
        
        diss_rwa = diss_eff[rwa]
        @show diss_rwa
        
        γ_rwa = γ_eff[rwa]
        @show γ_rwa
    end

    # -Δ, -Δ:
    begin
        @test isequal(simplify(Ω_rwa[3] - -2*Δ), 0)
        @test issetequal(diss_rwa[3], (a*σp, a*σp))

        @test issetequal(
            simplify.(expand.(γ_rwa[3].exponents)),
            simplify.(expand.([2*Δ^2, 4*Δ^2]))
        )

        @test isequal(
            simplify(γ_rwa[3].prefacs[1] - im*g^2/(2*Δ)),
            0
        )

        @test isequal(
            simplify(γ_rwa[3].prefacs[2] - -im*g^2/(2*Δ)),
            0
        )

    end

    # Δ, Δ:
    begin
        @test isequal(simplify(Ω_rwa[6] - 2*Δ), 0)
        @test issetequal(diss_rwa[6], (a'*σm, a'*σm))

        @test issetequal(
            simplify.(expand.(γ_rwa[6].exponents)),
            simplify.(expand.([2*Δ^2, 4*Δ^2]))
        )

        @test isequal(
            simplify(γ_rwa[6].prefacs[1] - -im*g^2/(2*Δ)),
            0
        )

        @test isequal(
            simplify(γ_rwa[6].prefacs[2] - im*g^2/(2*Δ)),
            0
        )
    end
    
end