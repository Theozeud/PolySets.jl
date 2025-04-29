function Base.:-(ps::PolySet{TP}) where {TP}
    qs = similar(ps)
    minus!(qs,ps)
    qs 
end

function Base.:+(ps::PolySet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    add!(qs, ps, c)
    qs 
end

function Base.:-(ps::PolySet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    minus!(qs,ps,c)
    qs 
end

function Base.:*(ps::PolySet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    mul!(qs,ps,c)
    qs 
end

function Base.:/(ps::PolySet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    div!(qs,ps,c)
    qs 
end

Base.:+(c::TC, ps::PolySet{TP}) where {TP, TC<:Number} = ps+c
Base.:-(c::TC, ps::PolySet{TP}) where {TP, TC<:Number} = ps+c
Base.:*(c::TC, ps::PolySet{TP}) where {TP, TC<:Number} = ps*c


function minus!(qs::PolySet, ps::PolySet)
    @. qs.coeffs = -ps.coeffs
    qs 
end

function add!(qs::PolySet, ps::PolySet, c::TC) where {TC<:Number}
    qs.coeffs .= ps.coeffs
    @views vqs = qs.coeffs[:,1]
    vqs .+= c 
    qs 
end


function minus!(qs::PolySet, ps::PolySet, c::TC) where {TC<:Number}
    qs.coeffs .= ps.coeffs
    @views vqs = qs.coeffs[:,1]
    vqs .-= c 
    qs 
end

function LinearAlgebra.mul!(qs::PolySet, ps::PolySet, c::TC) where {TC<:Number}
    @. qs.coeffs = ps.coeffs * c
    qs 
end

function div!(qs::PolySet, ps::PolySet, c::TC) where {TC<:Number}
    @. qs.coeffs = ps.coeffs / c
    qs 
end



#=
"""
    mul(ps::PolySet{TP}, qs::PolySet{TQ}) where {TP, TQ}

Returns a new `PolySet` containing the product of every polynomial in `ps` with every polynomial in `qs`.

# Arguments
- `ps::PolySet{TP}`: A set of `np` polynomials, each of degree `degp`.
- `qs::PolySet{TQ}`: A set of `nq` polynomials, each of degree `degq`.

# Returns
- `PolySet{NewT}`: A new `PolySet` of `np × nq` polynomials, each of degree at most `degp + degq`, 
   where `NewT` is the promoted type of `TP` and `TQ`.

# Details
This function constructs a temporary matrix `M` used to perform the polynomial multiplication via matrix multiplication. 
For each polynomial in `qs`, it fills `M` with shifted copies of its coefficients along diagonals, then computes all products with `ps` in one matrix-matrix multiplication.
"""
function Base.:*(ps::PolySet{TP}, qs::PolySet{TQ}) where {TP,TQ}
    (np, degp) = size(ps)
    (nq, degq) = size(qs)
    NewT = promote_type(TP,TQ)
    result = allocate_PolySet(NewT, np*nq, degp+degq-2)
    M = zeros(NewT,degp,degp+degq-1)
     @inbounds for i ∈ 1:nq
        @views coeffs = qs.coeffs[i,:]
        fill_upper_diagonals!(M, coeffs)
        @views vresult = result.coeffs[(i-1)*np+1:i*np,:]
        mul!(vresult, ps.coeffs, M)
    end
    result 
end

pairwiseproduct(ps::PolySet) = mul(ps,ps)



function factorize(ps::PolySet{TP}, coeffs::Vector{T}) where {TP,T}
    NewT = promote_type(TP,T)
    N = findlast(!iszero, coeffs)
    @assert !isnothing(N) "You can not factorize a polynomial by the null polynomial."
    degP = maxdeg(ps) - N +1
    @assert degP ≥ 0 "You can not factorise a polynomial by a polynomial with a strictly higher degree."
    P = allocate_PolySet(NewT, npolys(ps),degP+1)
    M = zeros(NewT, maxdeg(ps)+1 , maxdeg(ps)+1)
    @views vcoeffs = coeffs[1:N]
    fill_upper_diagonals2!(M, vcoeffs)
    @show size(P.coeffs)
    ps.coeffs / M

end
=#