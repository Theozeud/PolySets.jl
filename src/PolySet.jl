module PolySet

    using LinearAlgebra
    using SparseArrays
    using UnPack

    include("structure.jl")
    include("evaluate.jl")
    include("integrate.jl")
    include("derivate.jl")
end
