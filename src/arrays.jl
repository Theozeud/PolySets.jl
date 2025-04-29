"""
    sparse(ps::PolysSet)

Convert the coefficients matrix of ps into a sparse matrix.
"""
function SparseArrays.sparse(ps::PolysSet)
    sparseps = sparse(ps.coeffs)
    PolysSet(sparseps)
end