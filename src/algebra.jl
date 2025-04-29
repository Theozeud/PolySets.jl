function Base.:-(ps::PolysSet{TP}) where {TP}
    qs = similar(ps)
    minus!(qs,ps)
    qs 
end

function Base.:+(ps::PolysSet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    add!(qs, ps, c)
    qs 
end

function Base.:-(ps::PolysSet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    minus!(qs,ps,c)
    qs 
end

function Base.:*(ps::PolysSet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    mul!(qs,ps,c)
    qs 
end

function Base.:/(ps::PolysSet{TP}, c::TC) where {TP,TC<:Number}
    NewT = promote_type(TP,TC)
    qs = similar(ps, NewT)
    div!(qs,ps,c)
    qs 
end

Base.:+(c::TC, ps::PolysSet{TP}) where {TP, TC<:Number} = ps+c
Base.:-(c::TC, ps::PolysSet{TP}) where {TP, TC<:Number} = ps+c
Base.:*(c::TC, ps::PolysSet{TP}) where {TP, TC<:Number} = ps*c


function minus!(qs::PolysSet, ps::PolysSet)
    @. qs.coeffs = -ps.coeffs
    qs 
end

function add!(qs::PolysSet, ps::PolysSet, c::TC) where {TC<:Number}
    qs.coeffs .= ps.coeffs
    @views vqs = qs.coeffs[:,1]
    vqs .+= c 
    qs 
end


function minus!(qs::PolysSet, ps::PolysSet, c::TC) where {TC<:Number}
    qs.coeffs .= ps.coeffs
    @views vqs = qs.coeffs[:,1]
    vqs .-= c 
    qs 
end

function LinearAlgebra.mul!(qs::PolysSet, ps::PolysSet, c::TC) where {TC<:Number}
    @. qs.coeffs = ps.coeffs * c
    qs 
end

function div!(qs::PolysSet, ps::PolysSet, c::TC) where {TC<:Number}
    @. qs.coeffs = ps.coeffs / c
    qs 
end



#=
"""
    mul(ps::PolysSet{TP}, qs::PolysSet{TQ}) where {TP, TQ}

Returns a new `PolysSet` containing the product of every polynomial in `ps` with every polynomial in `qs`.

# Arguments
- `ps::PolysSet{TP}`: A set of `np` polynomials, each of degree `degp`.
- `qs::PolysSet{TQ}`: A set of `nq` polynomials, each of degree `degq`.

# Returns
- `PolysSet{NewT}`: A new `PolysSet` of `np × nq` polynomials, each of degree at most `degp + degq`, 
   where `NewT` is the promoted type of `TP` and `TQ`.

# Details
This function constructs a temporary matrix `M` used to perform the polynomial multiplication via matrix multiplication. 
For each polynomial in `qs`, it fills `M` with shifted copies of its coefficients along diagonals, then computes all products with `ps` in one matrix-matrix multiplication.
"""
function Base.:*(ps::PolysSet{TP}, qs::PolysSet{TQ}) where {TP,TQ}
    (np, degp) = size(ps)
    (nq, degq) = size(qs)
    NewT = promote_type(TP,TQ)
    result = allocate_PolysSet(NewT, np*nq, degp+degq-2)
    M = zeros(NewT,degp,degp+degq-1)
     @inbounds for i ∈ 1:nq
        @views coeffs = qs.coeffs[i,:]
        fill_upper_diagonals!(M, coeffs)
        @views vresult = result.coeffs[(i-1)*np+1:i*np,:]
        mul!(vresult, ps.coeffs, M)
    end
    result 
end

pairwiseproduct(ps::PolysSet) = mul(ps,ps)



function factorize(ps::PolysSet{TP}, coeffs::Vector{T}) where {TP,T}
    NewT = promote_type(TP,T)
    N = findlast(!iszero, coeffs)
    @assert !isnothing(N) "You can not factorize a polynomial by the null polynomial."
    degP = maxdeg(ps) - N +1
    @assert degP ≥ 0 "You can not factorise a polynomial by a polynomial with a strictly higher degree."
    P = allocate_PolysSet(NewT, npolys(ps),degP+1)
    M = zeros(NewT, maxdeg(ps)+1 , maxdeg(ps)+1)
    @views vcoeffs = coeffs[1:N]
    fill_upper_diagonals2!(M, vcoeffs)
    @show size(P.coeffs)
    ps.coeffs / M

end
=#