module PolySet

    using LinearAlgebra
    using SparseArrays
    using UnPack

    export  PolysSet, 
            allocate_PolysSet, npolys, maxdeg,
            evaluate, evaluate!,
            integrate, integrate!,
            derivate, derivate!

    include("structure.jl")
    include("evaluate.jl")
    include("integrate.jl")
    include("derivate.jl")

    export Monomials, Legendre, IntLegendre
    include("basis.jl")
end
