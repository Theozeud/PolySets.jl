module PolySet

    using LinearAlgebra
    using SparseArrays
    using UnPack
    using Plots

    export  PolysSet, 
            allocate_PolysSet, npolys, maxdeg,
            evaluate, evaluate!,
            integrate, integrate!,
            derivate, derivate!

    include("structure.jl")
    include("evaluate.jl")
    include("integrate.jl")
    include("derivate.jl")

    export minus!, add!, mul!, div!
    include("algebra.jl")

    export Monomials, Legendre, IntLegendre
    include("basis.jl")

    export sparse
    include("arrays.jl")

    export plot, plot!
    include("plot.jl")
end
