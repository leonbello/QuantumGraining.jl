using Richardson
using QuantumCumulants
using SymbolicUtils
using Symbolics

f(x) = sin(x)/x

extrapolate(1.0, rtol=1e-10) do x
    @show x
    sin(x)/x
end

extrapolate(1.0, x0=Inf) do x
    @show x
    (x^2 + 3x - 2) / (x^2 + 5)
end

extrapolate(1.0, x0=0) do x
    @show x
    1/x - 1/x
end

@cnumbers a x
f(x) = a*sin(x)/x + 1
extrapolate(1.0, x0=0) do x
    @show x
    f(x)
end