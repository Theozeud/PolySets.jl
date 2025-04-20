"""
    PolysSet{T, typeData<:AbstractMatrix}

A structure for efficient, vectorized manipulation of multiple polynomials.

# Fields
- `coeffs::typeData`: a matrix of polynomial coefficients. Each row corresponds to a single polynomial, with coefficients ordered from lowest to highest degree.

# Constructor
- `PolysSet(polys::AbstractMatrix)`: constructs a `PolysSet` from a matrix `polys`, where each row is interpreted as a polynomial. The type parameter `T` is inferred as the element type of the matrix, and `typeData` is the concrete matrix type (e.g., `Matrix{Float64}`, `SMatrix`, etc.).

# Examples
```julia
julia> A = [1 2 3; 0 1 0]
2×3 Matrix{Int64}:
 1  2  3
 0  1  0

julia> ps = PolysSet(A)
PolysSet{Int64, Matrix{Int64}}([1 2 3; 0 1 0])

``` 
"""
struct PolysSet{T, typeData<:AbstractMatrix}
    coeffs::typeData
    function PolysSet(polys::AbstractMatrix)
        new{eltype(polys), typeof(polys)}(polys)
    end
end



# Type element
Base.eltype(p::PolysSet) = eltype(p.coeffs)

# Indexing
Base.getindex(p::PolysSet, i::Int, j::Int) = p.coeffs[i, j]
Base.getindex(p::PolysSet, i::Int) = p.coeffs[i, :]
Base.setindex!(p::PolysSet, v, i::Int, j::Int) = (p.coeffs[i, j] = v)
Base.setindex!(p::PolysSet, v::AbstractVector, i::Int) = (p.coeffs[i, :] = v)
Base.first(ps::PolysSet) = ps[1, :]
Base.last(ps::PolysSet) = ps[end, :]
Base.firstindex(ps::PolysSet) = firstindex(ps.coeffs)
Base.lastindex(ps::PolysSet) = lastindex(ps.coeffs)


# Size
Base.size(p::PolysSet) = size(p.coeffs)
Base.size(p::PolysSet, dim::Int) = size(p.coeffs, dim)
Base.length(ps::PolysSet) = size(ps, 1)
npolys(ps::PolysSet) = size(ps.coeffs, 1)
maxdeg(ps::PolysSet) = size(ps.coeffs, 2) - 1

# Iteration
Base.iterate(p::PolysSet, state...) = iterate(eachrow(p.coeffs), state...)
Base.iterate(p::PolysSet) = iterate(p.coeffs)
Base.iterate(p::PolysSet, state) = iterate(p.coeffs, state)

# Display
Base.show(io::IO, p::PolysSet) = show(io, p.coeffs)

# Copy
Base.copy(p::PolysSet) = PolysSet(copy(p.coeffs))

# Similar and zeros
Base.similar(p::PolysSet) = PolysSet(similar(p.coeffs))
Base.similar(p::PolysSet, ::Type{T}) where T = PolysSet(similar(p.coeffs, T))
Base.similar(p::PolysSet, dims::Dims) = PolysSet(similar(p.coeffs, dims))
Base.similar(p::PolysSet, ::Type{T}, dims::Dims) where T = PolysSet(similar(p.coeffs, T, dims))

Base.zeros(::Type{PolysSet{T}}, dims::Dims) where {T} = PolysSet(zeros(T, dims))
Base.zeros(p::PolysSet) = PolysSet(zeros(eltype(p), size(p)))

# Convert and promote
Base.convert(::Type{PolysSet{T}}, A::AbstractMatrix{T}) where {T} = PolysSet{T, typeof(A)}(A)
Base.promote_rule(::Type{PolysSet{T}}, ::Type{PolysSet{S}}) where {T, S} = PolysSet{promote_type(T, S)}

# Equality
Base.:(==)(p1::PolysSet, p2::PolysSet) = p1.coeffs == p2.coeffs
Base.:(!=)(p1::PolysSet, p2::PolysSet) = !(p1 == p2)



"""
    allocate_PolysSet(T::DataType, nbpoly::Int, degmax::Int)

Create a`PolysSet` of `nbpoly` polynomials, each of degree at most `degmax` and whose 
coefficient matrix is filled with zeros.

# Arguments
- `T::DataType`: the type of the polynomial coefficients.
- `nbpoly::Int`: the number of polynomials to allocate.
- `degmax::Int`: the maximum degree of each polynomial.

# Returns
- A `PolysSet{T}` initialized with zero coefficients.
"""
allocate_PolysSet(T::DataType, nbpoly::Int, degmax::Int) = PolysSet(zeros(T,nbpoly, degmax+1))




function SparseArrays.sparse(ps::PolysSet)
    sparseps = sparse(ps.coeffs)
    PolysSet(sparseps)
end





function Base.show(io::IO, ::MIME"text/plain", ps::PolysSet)
    n, degmax_plus1 = size(ps.coeffs)
    println(io, "PolysSet with $n polynomials (degmax = $(degmax_plus1 - 1))")
    show(io, MIME"text/plain"(), ps.coeffs)
end