using Revise
using Symbolics
using Test
using QuantumCumulants
using QuantumGraining

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

# Calculate the TCG dissipators up to the second order
Iγ_list  = []
ω_list  = []
J_list  = []

for i in 1:2  
    perm_h, perm_g, perm_Ω = repeated_combinations(hvec, gvec, Ω, i)

    for (ω, g, h) in zip(perm_Ω, perm_g, perm_h) 
        for l in 1:(i-1)
            Iγ = (contraction_coeff(l, i-l, ω) - contraction_coeff(i-l, l, -reverse(ω)))
            Iγ = merge_duplicate_exponents(Iγ)
            J = prod(h[1:l]) # Jump operators
            L = prod(h[l+1:i])

            if !isequal(expand(sum(Iγ.prefacs.*Iγ.prefacs)),0)
                push!(Iγ_list, prod(g)*Iγ)
                push!(ω_list, sum(ω))
                push!(J_list, (J, L))
            end
        end
    end
end


@show J_list


# ((a*σm), (a*σm)): 2ωa + 2ωc

TestExpASigmaMinus_ASigmaMinus = [4((ωa + ωc)^2), 2((ωa + ωc)^2)]
TestPrefacASigmaMinus_ASigmaMinus = [[-(g^2//2)*1/(ωa + ωc)],
[(g^2//2)*1/(ωa + ωc)]
]

begin
    @test isequal(J_list[1], ((a*σm), (a*σm)))
    @test isequal(ω_list[1], 2*ωa + 2*ωc)

    @test isequal(length(Iγ_list[1].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[1].prefacs[i].*Iγ_list[1].polys[i][j]
                -
                TestPrefacASigmaMinus_ASigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[1].exponents[i], TestExpASigmaMinus_ASigmaMinus[i])
    end
end


# ((a*σp), (a*σm)): 2ωc

TestExpASigmaPlus_ASigmaMinus = [(ωa + ωc)^2 + (-ωa + ωc)^2 + 2(ωa + ωc)*(-ωa + ωc), (ωa + ωc)^2 + (-ωa + ωc)^2]
TestPrefacASigmaPlus_ASigmaMinus = [[(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))],
[-(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[2], ((a*σp), (a*σm)))
    @test isequal(ω_list[2], 2*ωc)

    @test isequal(length(Iγ_list[2].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[2].prefacs[i].*Iγ_list[2].polys[i][j]
                -
                TestPrefacASigmaPlus_ASigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[2].exponents[i], TestExpASigmaPlus_ASigmaMinus[i])
    end
end


# ((a'*σm), (a*σm)): 2ωa

TestExpAdSigmaMinus_ASigmaMinus = [(ωa + ωc)^2 + (ωa - ωc)^2 + 2(ωa + ωc)*(ωa - ωc), (ωa + ωc)^2 + (ωa - ωc)^2]
TestPrefacAdSigmaMinus_ASigmaMinus = [[-(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))],
[(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[3], ((a'*σm), (a*σm)))
    @test isequal(ω_list[3], 2*ωa)

    @test isequal(length(Iγ_list[3].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[3].prefacs[i].*Iγ_list[3].polys[i][j]
                -
                TestPrefacAdSigmaMinus_ASigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[3].exponents[i], TestExpAdSigmaMinus_ASigmaMinus[i])
    end
end


# ((a'*σp), (a'*σp)): - 2ωa - 2ωc

TestExpAdSigmaPlus_AdSigmaPlus = [4((-ωa - ωc)^2), 2((-ωa - ωc)^2)]
TestPrefacAdSigmaPlus_AdSigmaPlus = [[(g^2//2)*1/(ωa + ωc)],
[-(g^2//2)*1/(ωa + ωc)]
]

begin
    @test isequal(J_list[4], ((a'*σp), (a'*σp)))
    @test isequal(ω_list[4], - 2*ωa - 2*ωc)

    @test isequal(length(Iγ_list[4].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[4].prefacs[i].*Iγ_list[4].polys[i][j]
                -
                TestPrefacAdSigmaPlus_AdSigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[4].exponents[i], TestExpAdSigmaPlus_AdSigmaPlus[i])
    end
end


# ((a*σp), (a'*σp)): -2ωa

TestExpASigmaPlus_AdSigmaPlus = [(-ωa - ωc)^2 + (-ωa + ωc)^2 + 2(-ωa + ωc)*(-ωa - ωc), (-ωa - ωc)^2 + (-ωa + ωc)^2]
TestPrefacASigmaPlus_AdSigmaPlus = [[(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))],
[-(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[5], ((a*σp), (a'*σp)))
    @test isequal(ω_list[5], - 2*ωa)

    @test isequal(length(Iγ_list[5].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[5].prefacs[i].*Iγ_list[5].polys[i][j]
                -
                TestPrefacASigmaPlus_AdSigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[5].exponents[i], TestExpASigmaPlus_AdSigmaPlus[i])
    end
end


# ((a'*σge), (a'*σeg)): -2ωc

TestExpAdSigmaMinus_AdSigmaPlus = [(ωa - ωc)^2 + (-ωa - ωc)^2 + 2(ωa - ωc)*(-ωa - ωc), (ωa - ωc)^2 + (-ωa - ωc)^2]
TestPrefacAdSigmaMinus_AdSigmaPlus = [[-(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))],
[(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[6], ((a'*σm), (a'*σp)))
    @test isequal(ω_list[6], - 2*ωc)

    @test isequal(length(Iγ_list[6].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[6].prefacs[i].*Iγ_list[6].polys[i][j]
                -
                TestPrefacAdSigmaMinus_AdSigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[6].exponents[i], TestExpAdSigmaMinus_AdSigmaPlus[i])
    end
end


# ((a*σm), (a*σp)): 2ωc

TestExpASigmaMinus_ASigmaPlus = [(ωa + ωc)^2 + (-ωa + ωc)^2 + 2(ωa + ωc)*(-ωa + ωc), (ωa + ωc)^2 + (-ωa + ωc)^2]
TestPrefacASigmaMinus_ASigmaPlus = [[(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))],
[-(g^2//2)*ωc/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[7], ((a*σm), (a*σp)))
    @test isequal(ω_list[7], 2*ωc)

    @test isequal(length(Iγ_list[7].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[7].prefacs[i].*Iγ_list[7].polys[i][j]
                -
                TestPrefacASigmaMinus_ASigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[7].exponents[i], TestExpASigmaMinus_ASigmaPlus[i])
    end
end


# ((a'*σp), (a*σp)): -2ωa

TestExpAdSigmaPlus_ASigmaPlus = [(-ωa - ωc)^2 + (-ωa + ωc)^2 + 2(-ωa + ωc)*(-ωa - ωc), (-ωa - ωc)^2 + (-ωa + ωc)^2]
TestPrefacAdSigmaPlus_ASigmaPlus = [[(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))],
[-(g^2//2)*ωa/((ωa + ωc)*(ωa - ωc))]
]

begin
    @test isequal(J_list[8], ((a'*σp), (a*σp)))
    @test isequal(ω_list[8], -2*ωa)

    @test isequal(length(Iγ_list[8].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[8].prefacs[i].*Iγ_list[8].polys[i][j]
                -
                TestPrefacAdSigmaPlus_ASigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[8].exponents[i], TestExpAdSigmaPlus_ASigmaPlus[i])
    end
end


# ((a*σp), (a*σp)): −2ωa+2ωc

TestExpASigmaPlus_ASigmaPlus = [4((-ωa + ωc)^2), 2((-ωa + ωc)^2)]
TestPrefacASigmaPlus_ASigmaPlus = [[(g^2//2)*1/(ωa - ωc)],
[-(g^2//2)*1/(ωa - ωc)]
]

begin
    @test isequal(J_list[9], ((a*σp), (a*σp)))
    @test isequal(ω_list[9], 2*ωc − 2*ωa)

    @test isequal(length(Iγ_list[9].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[9].prefacs[i].*Iγ_list[9].polys[i][j]
                -
                TestPrefacASigmaPlus_ASigmaPlus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[9].exponents[i], TestExpASigmaPlus_ASigmaPlus[i])
    end
end


# ((a*σm), (a'*σm)): 2ωa

TestExpASigmaMinus_AdSigmaMinus = [(ωa + ωc)^2 + (ωa - ωc)^2 + 2(ωa + ωc)*(ωa - ωc), (ωa + ωc)^2 + (ωa - ωc)^2]
TestPrefacASigmaMinus_AdSigmaMinus = [[-(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))],
[(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))]
]

begin
    @test isequal(J_list[10], ((a*σm), (a'*σm)))
    @test isequal(ω_list[10], 2*ωa)

    @test isequal(length(Iγ_list[10].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[10].prefacs[i].*Iγ_list[10].polys[i][j]
                -
                TestPrefacASigmaMinus_AdSigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[10].exponents[i], TestExpASigmaMinus_AdSigmaMinus[i])
    end
end


# ((a'*σp), (a'*σm)): -2ωc

TestExpAdSigmaPlus_AdSigmaMinus = [(ωa - ωc)^2 + (-ωa - ωc)^2 + 2(ωa - ωc)*(-ωa - ωc), (ωa - ωc)^2 + (-ωa - ωc)^2]
TestPrefacAdSigmaPlus_AdSigmaMinus = [[-(g^2//2)*ωc/((ωa - ωc)*(ωa + ωc))],
[(g^2//2)*ωc/((ωa - ωc)*(ωa + ωc))]
]


begin
    @test isequal(J_list[11], ((a'*σp), (a'*σm)))
    @test isequal(ω_list[11], -2*ωc)

    @test isequal(length(Iγ_list[11].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[11].prefacs[i].*Iγ_list[11].polys[i][j]
                -
                TestPrefacAdSigmaPlus_AdSigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[11].exponents[i], TestExpAdSigmaPlus_AdSigmaMinus[i])
    end
end


# ((a'*σm), (a'*σm)): 2ωa−2ωc

TestExpAdSigmaPlus_AdSigmaMinus = [4((ωa - ωc)^2), 2((ωa - ωc)^2)]
TestPrefacAdSigmaPlus_AdSigmaMinus = [[-(g^2//2)*1/(ωa - ωc)],
[(g^2//2)*1/(ωa - ωc)]
]


begin
    @test isequal(J_list[12], ((a'*σm), (a'*σm)))
    @test isequal(ω_list[12], 2*ωa−2*ωc)

    @test isequal(length(Iγ_list[12].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                Iγ_list[12].prefacs[i].*Iγ_list[12].polys[i][j]
                -
                TestPrefacAdSigmaPlus_AdSigmaMinus[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(Iγ_list[12].exponents[i], TestExpAdSigmaPlus_AdSigmaMinus[i])
    end
end


@test true