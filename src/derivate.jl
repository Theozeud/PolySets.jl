"""
    derivate(ps::PolySet) -> PolySet

Returns a new `PolySet` containing the derivatives of the polynomials in `ps`.

Differentiation is performed in a fully vectorized way, and the result is returned as a new `PolySet`.
"""
function derivate(ps::PolySet)
    ips = similar(ps, (npolys(ps),maxdeg(ps)))
    derivate!(ips, ps)
    ips
end


"""
    derivate!(ips::PolySet, ps::PolySet)

Computes the derivative of each polynomial in `ps` and stores the result in `ips` (in-place).

The matrix `ips.coeffs` must have the same number of rows as `ps`, and at least `maxdeg(ps)` columns
(since the result of a derivative has one fewer term). The operation is fully vectorized across rows.
"""
function derivate!(ips::PolySet, ps::PolySet)
    fill!(ips.coeffs, 0)
    @views vips = ips.coeffs[:,1:maxdeg(p)]
    vips .= ps[:,2:end] 
    x = 1:maxdeg(ps)
    vips = vips .* x'
end