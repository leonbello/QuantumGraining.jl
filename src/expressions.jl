"""
    @definemodes

Helper macro to define `n` up-modes and `m` down-modes. Generates two vectors of symbols.
        
# Examples
    @definemodes μ 3 ν 5    
    > μ = [μ1 μ2 μ3]
    > ν = [ν1 ν2 ν3 ν4 ν5]
    > typeof(μ) = Array{QuantumCumulants.cnumber}
"""
macro definemodes(uname, unum)
    esc_uname = esc(uname)
    esc_unum = esc(unum)
    uexp, u_arr = _definemodes(uname, eval(unum))
    return quote
        @cnumbers τ
        $uexp
        $esc_uname = eval.($(u_arr))
    end
end

"""
    _definemodes(name, n)

Returns an expression of the form @cnumbers name1 name2...
"""
function _definemodes(name::Base.Symbol, n::Int, fmt=:QuantumCumulants)
    new_ex = :(@cnumbers)
    sym_arr = []
    for i in 1:n
        ex = Symbol(name, i)
        sym_arr = push!(sym_arr, ex)
        push!(new_ex.args, ex)
    end
    return (new_ex, sym_arr)
end


# macro definemodes(n, m) 

#     umodes = eval(:([Symbol(:μ, i) for i in range(1, $n)]))
#     ex1 = :(@cnumbers)
#     for mode in umodes
#         push!(ex1.args, mode)
#     end                                                                     # @cnumbers μ1 μ2 ... μn

#     dmodes = [Symbol(:ν, i) for i in range(1, eval(:($m)) )]
#     ex2 = :(@cnumbers)
#     for mode in dmodes
#         push!(ex2.args, mode)
#     end                                                                     # @cnumbers ν1 ν2 ... νm
    
#     ex3 = quote
#         μ = [eval(Symbol(:μ, i)) for i in range(1, $n)]
#         ν = [eval(Symbol(:ν, i)) for i in range(1, $m)]
#     end                                                                     # generates a list of the mode variables

#     ex4 = quote
#         @cnumbers τ  
#     end
#     return esc(:($ex1; $ex2; $ex3; $ex4))
# end

# macro definemodes(sym1, sym2, diagram)
#     #sym1 = esc(sym1)
#     #sym2 = esc(sym2)                                                        # take the names the user gives
#     n, m = eval(:(_maxmodes($diagram)))

#     umodes = [Symbol(sym1, i) for i in range(1, n)]
#     #umodes = eval( :([Symbol(sym1, i) for i in range(1, $n)]) )
#     ex1 = :(@cnumbers)
#     for mode in umodes
#         push!(ex1.args, mode)
#     end                                                                     # @cnumbers μ1 μ2 ... μn

#     dmodes = [Symbol(sym2, i) for i in range(1, m)]
#     ex2 = :(@cnumbers)
#     for mode in dmodes
#         push!(ex2.args, mode)
#     end                                                                     # @cnumbers ν1 ν2 ... νm
    
#     μ = Symbol(sym1, "_arr")
#     ν = Symbol(sym2, "_arr")
#     ex3 = quote
#         $(μ) = [eval(Symbol($sym1, i)) for i in range(1, $n)]
#         $(ν) = [eval(Symbol($sym2, i)) for i in range(1, $m)]
#     end                                                                     # generates a list of the mode variables

#     ex4 = quote
#         @cnumbers τ  
#     end
#     return esc(:($ex1; $ex2; $ex3; $ex4))
# end