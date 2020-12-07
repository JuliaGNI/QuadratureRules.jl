"""
Clenshaw-Curtis quadrature.
"""
function ClenshawCurtisQuadrature(::Type{T}, s::Integer; IT=BigFloat) where {T}
    if s == 1
        throw(ErrorException("Clenshaw-Curtins quadrature is not defined for one stage."))
    end

    mymapreduce(f, op, itr) = isempty(itr) ? zero(IT) : mapreduce(f, op, itr)

    b(j,n) = fld(n,2) == j == cld(n,2) ? IT(1) : IT(2)
    c(k,n) = mod(k,n) == 0 ? IT(1) : IT(2)
    ϑ(k,n) = @big 2π * k / n

    cctrm(k,j,n) = @big b(j,n) / ( 4 * j^2 - 1 ) * cos(j * ϑ(k,n))
    ccsum(k,n) = c(k,n) / IT(2n) * ( IT(1) - mymapreduce(j -> cctrm(k,j,n), +, 1:div(n,2)) )

    x = chebyshev_points(IT, s, Val(2))
    w = [ccsum(i-1,s-1) for i in 1:s]

    QuadratureRule(s^2, x, w, T)
end

ClenshawCurtisQuadrature(s) = ClenshawCurtisQuadrature(Float64, s)
