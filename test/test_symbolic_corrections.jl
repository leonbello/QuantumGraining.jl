using Revise
using IterTools
using Symbolics
using Test
using Revise
using QuantumGraining

begin
    @variables ω1 ω2 ω3 ω4 ω5
    μ1 = [ω4]
    ν1 = [ω5, -ω5]

    μ2 = [ω1, ω2, ω3, -ω3]
    ν2 = Num[];

    ω = [(μ1, ν1), (μ2, ν2)]

    up_poles, down_poles = find_all_poles(ω);
    @show down_poles
    @test up_poles == [Int[], [2]]
    @test down_poles == [[2], []]
    @test count_poles(find_all_poles(ω)...) == 2

    corr = diagram_correction(ω)
    @show corr
end

begin
    @variables μ11 ν11 ν12 μ21 μ22 μ23 μ24

    μ1 = UVec([μ11])
    ν1 = DVec([ν11, ν12])
    μ2 = UVec([μ21, μ22, μ23, μ24])
    ν2 = DVec(Num[]);

    ω = [(μ1, ν1), (μ2, ν2)]

    up_poles = Int[]
    down_poles = [2]
    u, d = 1, 10
    calc_pole_corrections(μ1, ν1, up_poles, down_poles, u, d)
end