"""
    integrate(ps::PolySet{TPS}, a::Real, b::Real) where TPS

Evaluate the definite integral of each polynomial in the `PolySet` `ps` over the interval [`a`, `b`]. 

The result `y[i]` is the integral of the `i`-th polynomial from `a` to `b`.

# Arguments
- `ps`: A `PolySet` structure containing the coefficient matrix.
- `a`, `b`: Integration bounds.
"""

function integrate(ps::PolySet{TPS}, a::Real, b::Real)
    NewT = promote_type(TPS, typeof(a), typeof(b))
    y = zeros(NewT, npolys(ps))
    integrate!(y, ps, a, b)
    y
end


"""
    integrate!(y::AbstractVector, ps::PolySet{TPS}, a::Real, b::Real) where TPS

Evaluate the definite integral of each polynomial in the `PolySet` `ps` over the interval [`a`, `b`],
and store the result in the output vector `y`. 

The result `y[i]` is the integral of the `i`-th polynomial from `a` to `b`.

This version allocates a temporary vector internally.

# Arguments
- `y`: Output vector, must have length equal to the number of polynomials in `ps`.
- `ps`: A `PolySet` structure containing the coefficient matrix.
- `a`, `b`: Integration bounds.
"""
function integrate!(y::AbstractVecOrMat, ps::PolySet{TPS}, a::Real, b::Real) where TPS
    @unpack coeffs = ps
    (_,m) = size(ps)
    newT = promote_type(TPS, typeof(a), typeof(b))
    if a != -b
        cache = zeros(newT,m)
        integrate!(y, ps, a, b, cache)
    else     
        s = div(m+1,2)
        cache = zeros(newT,s)
        integrate!(y, ps, a, b, cache)
    end
end


"""
    integrate!(y::AbstractVector, ps::PolySet{TPS}, a::Real, b::Real, cache::AbstractVector) where TPS

Same as `integrate!(y, ps, a, b)` but avoids allocation by reusing the provided `cache` vector.

# Arguments
- `y`: Output vector, must have length equal to the number of polynomials in `ps`.
- `ps`: A `PolySet` structure containing the coefficient matrix.
- `a`, `b`: Integration bounds.
- `cache`: Pre-allocated vector of length equal to the polynomial degree, used to store powers.

"""
function integrate!(y::AbstractVecOrMat, ps::PolySet{TPS}, a::Real, b::Real, cache::AbstractVector) where TPS
    @unpack coeffs = ps
    (_,m) = size(ps)
    newT = promote_type(TPS, typeof(a), typeof(b))
    if a != -b
        apow = one(newT)
        bpow = one(newT)
        @views vcache = cache[1:m]
        @inbounds for j in 1:m
            apow *= a
            bpow *= b
            vcache[j] = (bpow - apow) / j    
        end
        mul!(y,coeffs,vcache)
        return y
    else     
        s = div(m+1,2)
        @views vcache = cache[1:s]
        bpow = one(newT)
        @inbounds for j in 1:s
            bpow *= b
            vcache[j] = 2*bpow / (2*j -1)
            bpow *= b
        end
        @views coeffsv = coeffs[:,1:2:m]
        mul!(y,coeffsv,vcache)
        return y
    end
end


"""
    integrate!(ips::PolySet, ps::PolySet{TPS}, a::Real) where TPS

Computes the indefinite integral of each polynomial in the `PolySet` `ps`, evaluated to be zero at `x = a`, 
and stores the result in `ips`.

# Arguments
- `ips::PolySet`: A preallocated `PolySet` where the result will be stored. It must have the same number of 
polynomials as `ps` and one more degree.
- `ps::PolySet{TPS}`: A set of polynomials represented as a `PolySet` with scalar type `TPS`.
- `a::Real`: The point at which the antiderivative is set to zero. 

# Returns
- `ips::PolySet`: The modified `ips`, containing the integrated polynomials.

# Example
```julia
ps = PolySet([[1.0, 2.0], [0.0, 3.0]])             
ips = allocate_polyset(Float64, 2, 3)               
integrate!(ips, ps, 0.0)                                         
```
"""
function integrate!(ips::PolySet, ps::PolySet{TPS}, a::Real) where TPS
    x = 1 ./(1:(maxdeg(ps)+1))
    @views vips = ips.coeffs[:,2:maxdeg(ps)+2]
    mul!(vips, ps.coeffs, Diagonal(x))
    NewT = promote_type(TPS, typeof(a))
    y = zeros(NewT, npolys(ps))
    evaluate!(y, ips, a)
    @views vips = ips.coeffs[:,1]
    vips .-= y
    ips
end