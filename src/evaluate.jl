"""
    evaluate(ps::PolysSet, x::Number)

Evaluates a set of polynomials represented by `ps` at the scalar point `x`.

# Arguments
- `ps`: A `PolysSet` object storing polynomial coefficients. Each row corresponds to one polynomial; columns represent increasing powers of `x`.
- `x`: The scalar value at which all polynomials are evaluated.

# Notes
- Uses Horner's method for efficient polynomial evaluation.
"""
function evaluate(ps::PolysSet, x::Number)
    NewT = promote_type(typeof(x), eltype(ps))
    y = zeros(NewT, size(ps,1))
    evaluate!(y, ps, x)
end

(ps::PolysSet)(x::Number) = evaluate(ps, x)

"""
    evaluate(ps::PolysSet, x::Number)

Evaluates a set of polynomials represented by `ps` at the scalar point `x`.

# Arguments
- `ps`: A `PolysSet` object storing polynomial coefficients. Each row corresponds to one polynomial; columns represent increasing powers of `x`.
- `x`: A vector of scalar inputs at which the polynomials are evaluated.

# Notes
- Uses Horner's method for efficient polynomial evaluation.
"""
function evaluate(ps::PolysSet, x::AbstractVector)
    NewT = promote_type(eltype(x), eltype(ps))
    y = zeros(NewT, size(ps,1), length(x))
    evaluate!(y, ps, x)
end

(ps::PolysSet)(x::AbstractVector) = evaluate(ps, x)

"""
    evaluate!(y::AbstractVector, ps::PolysSet, x::Number)

Evaluates a set of polynomials represented by `ps` at the scalar point `x`, writing the result in-place into the vector `y`.

# Arguments
- `y`: Output vector (modified in-place), with length equal to the number of polynomials. After evaluation, `y[i]` contains the value of the i-th polynomial at `x`.
- `ps`: A `PolysSet` object storing polynomial coefficients. Each row corresponds to one polynomial; columns represent increasing powers of `x`.
- `x`: The scalar value at which all polynomials are evaluated.

# Notes
- Uses Horner's method for efficient polynomial evaluation.
"""
function evaluate!(y::AbstractVector, ps::PolysSet, x::Number)
    @unpack coeffs = ps
    @views y .= coeffs[:,end]
    @inbounds for i ∈ maxdeg(ps):-1:1
        @views coeffsi = coeffs[:,i]
        @. y = y * x + coeffsi
    end
    y
end

"""
    evaluate!(y::AbstractMatrix, ps::PolysSet, x::AbstractVector)

Evaluates a set of polynomials represented by `ps` at multiple scalar points given in `x`, writing the result in-place into the matrix `y`.

# Arguments
- `y`: Output matrix (modified in-place), of size `(n, m)` where `n` is the number of polynomials and `m == length(x)`. Column `j` of `y` contains the evaluations of all polynomials at `x[j]`.
- `ps`: A `PolysSet` object storing polynomial coefficients. Each row is a polynomial; columns represent increasing powers of `x`.
- `x`: A vector of scalar inputs at which the polynomials are evaluated.

# Notes
- Uses Horner's method for fast evaluation.
- Evaluation is vectorized over all input points.
"""
function evaluate!(y::AbstractMatrix, ps::PolysSet, x::AbstractVector)
    @unpack coeffs = ps
    fill!(y, 0)
    @views pend = coeffs[:,end]
    y .+= pend
    @inbounds for i ∈ maxdeg(ps):-1:1
        @views coeffsi = coeffs[:,i]
        @. y = y * x' + coeffsi
    end
    y
end
