#using QuantumGraining
using QuantumCumulants

module Tst
    include("../src/expressions.jl")

    function ex_test(n)
        @cnumbers a b c
        ex = n*a*b*c
        return ex
    end

    ex = ex_test(1)*ex_test(2)

    function ex_func(n)
        @cnumbers a
        func(n) = n*a
        return func
    end

    func = ex_func(3)
    func(3)
    
    @cnumbers ω1 ω2 ω3
    ex1 = ω1 + ω2
    ex2 = ω2 + ω3

    ex3 = ex1 - ex2

    # error when calling macro inside a function
    # solution one - don't use macros:
    function f1(x, n)
        ex, sym_arr = _definemodes(:x, n)
        eval(ex)
        print(x)
    end
    f1(x, 3)

    # solution two - figure out what's wrong with the macro
    # basically, a macro is just a way to automate pasting code -- it cannot depend on function arguments
    diagram3 = [(3,3), (2,1)]
    c3 = calculate_coeff(diagram3)

    function test_definemodes(diagram::Array{Tuple{Int, Int}})
        @definemodes diagram
    end
end
