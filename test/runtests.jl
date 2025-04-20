using PolySet
using Test

@testset "PolySet.jl" begin
    
    p = Legendre(10)
    evaluate(p, 2)
end
