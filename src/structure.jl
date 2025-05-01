"""
    PolySet{T, typeData<:AbstractMatrix}

A structure for efficient, vectorized manipulation of multiple polynomials.

# Fields
- `coeffs::typeData`: a matrix of polynomial coefficients. Each row corresponds to a single polynomial, with coefficients ordered from lowest to highest degree.

# Constructor
- `PolySet(polys::AbstractMatrix)`: constructs a `PolySet` from a matrix `polys`, where each row is interpreted as a polynomial. The type parameter `T` is inferred as the element type of the matrix, and `typeData` is the concrete matrix type (e.g., `Matrix{Float64}`, `SMatrix`, etc.).

# Examples
```julia
julia> A = [1 2 3; 0 1 0]
2×3 Matrix{Int64}:
 1  2  3
 0  1  0

julia> ps = PolySet(A)
PolySet{Int64, Matrix{Int64}}([1 2 3; 0 1 0])

``` 
"""
struct PolySet{T, typeData<:AbstractVecOrMat}
    coeffs::typeData
    function PolySet(polys::AbstractVecOrMat)
        new{eltype(polys), typeof(polys)}(polys)
    end
end


# Type element
Base.eltype(p::PolySet) = eltype(p.coeffs)

# Indexing
Base.getindex(p::PolySet, i::Int, j::Int) = p.coeffs[i, j]
Base.getindex(p::PolySet, I::UnitRange, j::Int) = p.coeffs[I, j]
Base.getindex(p::PolySet, i::Int, J::UnitRange) = p.coeffs[i, J]
Base.getindex(p::PolySet, I::UnitRange, J::UnitRange) = p.coeffs[I, J]
Base.getindex(p::PolySet, i::Int) = p.coeffs[i, :]
Base.getindex(p::PolySet, I::UnitRange) = p.coeffs[I, :]
Base.setindex!(p::PolySet, v, i::Int, j::Int) = (p.coeffs[i, j] = v)
Base.setindex!(p::PolySet, v::AbstractVector, i::Int) = (p.coeffs[i, :] = v)
Base.first(ps::PolySet) = ps[1, :]
Base.last(ps::PolySet) = ps[end, :]
Base.firstindex(ps::PolySet) = firstindex(ps.coeffs)
Base.lastindex(ps::PolySet) = lastindex(ps.coeffs)


# Size
Base.size(p::PolySet) = size(p.coeffs)
Base.size(p::PolySet, dim::Int) = size(p.coeffs, dim)
Base.length(ps::PolySet) = size(ps, 1)
npolys(ps::PolySet) = size(ps.coeffs, 1)
maxdeg(ps::PolySet) = size(ps.coeffs, 2) - 1

# Iteration
Base.iterate(p::PolySet, state...) = iterate(eachrow(p.coeffs), state...)
Base.iterate(p::PolySet) = iterate(p.coeffs)
Base.iterate(p::PolySet, state) = iterate(p.coeffs, state)

# Copy
Base.copy(p::PolySet) = PolySet(copy(p.coeffs))

# Similar and zeros
Base.similar(p::PolySet) = PolySet(similar(p.coeffs))
Base.similar(p::PolySet, ::Type{T}) where T = PolySet(similar(p.coeffs, T))
Base.similar(p::PolySet, dims::Dims) = PolySet(similar(p.coeffs, dims))
Base.similar(p::PolySet, ::Type{T}, dims::Dims) where T = PolySet(similar(p.coeffs, T, dims))

Base.zeros(::Type{PolySet{T}}, dims::Dims) where {T} = PolySet(zeros(T, dims))
Base.zero(p::PolySet) = PolySet(zeros(eltype(p), size(p)))

# Convert and promote
Base.convert(::Type{PolySet{T}}, A::AbstractMatrix{T}) where {T} = PolySet{T, typeof(A)}(A)
Base.promote_rule(::Type{PolySet{T}}, ::Type{PolySet{S}}) where {T, S} = PolySet{promote_type(T, S)}

# Equality
Base.:(==)(p1::PolySet, p2::PolySet) = p1.coeffs == p2.coeffs
Base.:(!=)(p1::PolySet, p2::PolySet) = !(p1 == p2)

# Display
function Base.show(io::IO, ::MIME"text/plain", ps::PolySet)
    n, degmax_plus1 = size(ps.coeffs)
    println(io, "PolySet with $n polynomials (degmax = $(degmax_plus1 - 1))")
    show(io, MIME"text/plain"(), ps.coeffs)
end


"""
    allocate_PolySet(T::DataType, nbpoly::Int, degmax::Int)

Create a`PolySet` of `nbpoly` polynomials, each of degree at most `degmax` and whose 
coefficient matrix is filled with zeros.

# Arguments
- `T::DataType`: the type of the polynomial coefficients.
- `nbpoly::Int`: the number of polynomials to allocate.
- `degmax::Int`: the maximum degree of each polynomial.

# Returns
- A `PolySet{T}` initialized with zero coefficients.
"""
allocate_PolySet(T::DataType, nbpoly::Int, degmax::Int) = PolySet(zeros(T,nbpoly, degmax+1))




