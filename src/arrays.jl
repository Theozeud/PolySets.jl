"""
    sparse(ps::PolySet)

Convert the coefficients matrix of ps into a sparse matrix.
"""
function SparseArrays.sparse(ps::PolySet)
    sparseps = sparse(ps.coeffs)
    PolySet(sparseps)
end