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

g_eff_2, Ω_eff_2 = effective_hamiltonian(hvec, gvec, Ω, 2; as_dict=true)

@show keys(g_eff_2)

# σee: 0

TestExpSigmaEE = [0, (ωa + ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + ωc)^2]
TestPrefacSigmaEE = [[(g^2//4)*(1/(ωa - ωc) + 1/(ωa + ωc))],
[(-(g^2//4)/(ωa + ωc))],
[(-(g^2//4)/(ωa - ωc))]
]

begin
    @test isequal(Ω_eff_2[σee], 0)

    @test isequal(length(g_eff_2[σee].exponents), 3)

    for i in 1:3
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[σee].prefacs[i].*g_eff_2[σee].polys[i][j]
                -
                TestPrefacSigmaEE[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[σee].exponents[i], TestExpSigmaEE[i])
    end
end


# a*a*σee: 2*ωc

TestExpAASigmaEE = [(4//1)*(ωc^2), (ωa + ωc)^2 + (-ωa + ωc)^2]
TestPrefacAASigmaEE = [[(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))],
[-(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))]
]

begin
    @test isequal(Ω_eff_2[a*a*σee], 2*ωc)

    @test isequal(length(g_eff_2[a*a*σee].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a*a*σee].prefacs[i].*g_eff_2[a*a*σee].polys[i][j]
                -
                TestPrefacAASigmaEE[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a*a*σee].exponents[i], TestExpAASigmaEE[i])
    end
end


# a'*a'σee: -2*ωc

TestExpAdAdSigmaEE = [(4//1)*(ωc^2), (ωa - ωc)^2 + (-ωa - ωc)^2]
TestPrefacAdAdSigmaEE = [[(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))],
[-(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))]
]

begin
    @test isequal(Ω_eff_2[a'*a'*σee], -2*ωc)

    @test isequal(length(g_eff_2[a'*a'*σee].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a'*a'*σee].prefacs[i].*g_eff_2[a'*a'*σee].polys[i][j]
                -
                TestPrefacAdAdSigmaEE[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a'*a'*σee].exponents[i], TestExpAdAdSigmaEE[i])
    end
end


# a'*a*σee: 0

TestExpAdASigmaEE = [0, (ωa + ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + ωc)^2]
TestPrefacAdASigmaEE = [[(g^2)*ωa/((ωa - ωc)*(ωa + ωc))],
[-(g^2//2)*1/(ωa + ωc)],
[-(g^2//2)*1/(ωa - ωc)]
]

begin
    @test isequal(Ω_eff_2[a'*a*σee], 0)

    @test isequal(length(g_eff_2[a'*a*σee].exponents), 3)

    for i in 1:3
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a'*a*σee].prefacs[i].*g_eff_2[a'*a*σee].polys[i][j]
                -
                TestPrefacAdASigmaEE[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a'*a*σee].exponents[i], TestExpAdASigmaEE[i])
    end
end


# a*a: 2*ωc

TestExpAA = [(4//1)*(ωc^2), (ωa + ωc)^2 + (-ωa + ωc)^2]
TestPrefacAA = [[-(g^2//4)*ωa/((ωa - ωc)*(ωa + ωc))],
[(g^2//4)*ωa/((ωa - ωc)*(ωa + ωc))]
]

begin
    @test isequal(Ω_eff_2[a*a], 2*ωc)

    @test isequal(length(g_eff_2[a*a].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a*a].prefacs[i].*g_eff_2[a*a].polys[i][j]
                -
                TestPrefacAA[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a*a].exponents[i], TestExpAA[i])
    end
end


# a'*a': -2*ωc

TestExpAdAd = [(4//1)*(ωc^2), (ωa - ωc)^2 + (-ωa - ωc)^2]
TestPrefacAdAd = [[-(g^2//4)*ωa/((ωa - ωc)*(ωa + ωc))],
[(g^2//4)*ωa/((ωa - ωc)*(ωa + ωc))]
]

begin
    @test isequal(Ω_eff_2[a'*a'], -2*ωc)

    @test isequal(length(g_eff_2[a'*a'].exponents), 2)

    for i in 1:2
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a'*a'].prefacs[i].*g_eff_2[a'*a'].polys[i][j]
                -
                TestPrefacAdAd[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a'*a'].exponents[i], TestExpAdAd[i])
    end
end


# a'*a: 0

TestExpAdA = [0, (ωa + ωc)^2 + (-ωa - ωc)^2, (ωa - ωc)^2 + (-ωa + ωc)^2]
TestPrefacAdA = [[-(g^2//2)*ωa/((ωa - ωc)*(ωa + ωc))],
[(g^2//4)*1/(ωa + ωc)],
[(g^2//4)*1/(ωa - ωc)]
]

begin
    @test isequal(Ω_eff_2[a'*a], 0)

    @test isequal(length(g_eff_2[a'*a].exponents), 3)

    for i in 1:3
        PrefacDiff = []
        for j in 1:1
            push!(PrefacDiff, simplify(
                g_eff_2[a'*a].prefacs[i].*g_eff_2[a'*a].polys[i][j]
                -
                TestPrefacAdA[i][j]
                ))
        end
        @test isequal(PrefacDiff, fill(0,1))
        @test isequal(g_eff_2[a'*a].exponents[i], TestExpAdA[i])
    end
end




# Test the 3rd-order Hamiltonian corrections with ωa=7/17, ωc=11/13, and g = 5/3

g_eff_3, Ω_eff_3 = effective_hamiltonian(hvec, gvec, Ω, 3; as_dict=true)

@show keys(g_eff_3)


# (a'*a'*a'*σp): −ωa−3ωc

TestExpAdAdAdSigmaP = [(-ωa - 3ωc)^2, (-ωa - ωc)^2 + (4//1)*(ωc^2), (ωa - ωc)^2 + 2((-ωa - ωc)^2)]
TestPrefacAdAdAdSigmaP = [[0.271847],
[-1.33092],
[1.05907]
]

begin
    @test isequal(Ω_eff_3[(a'*a'*a'*σp)], -ωa-3*ωc)

    @test isequal(length(g_eff_3[(a'*a'*a'*σp)].exponents), 3)

    for i in 1:3
        for j in 1:1
            diff = (substitute(g_eff_3[(a'*a'*a'*σp)].prefacs[i].*g_eff_3[(a'*a'*a'*σp)].polys[i][j], Dict(ωa => 7//17, ωc => 11//13, g => 5//3))
            -
            TestPrefacAdAdAdSigmaP[i][j])
            @test abs(diff/TestPrefacAdAdAdSigmaP[i][j]) < 0.0001
        end
        @test isequal(g_eff_3[(a'*a'*a'*σp)].exponents[i], TestExpAdAdAdSigmaP[i])
    end
end

@test true