
"""
Lobatto-Legendre points on the interval [-1,+1].

The points are computed in the arithmetic `IT` and converted to `T`.
"""
function lobatto_legendre_points(::Type{T}, s::Integer; IT=BigFloat) where {T}
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    D = Polynomials.derivative(Polynomial(IT[1, 0, -1])^(s-1), s-2)
    x = sort(_newton_roots(D, FastGaussQuadrature.gausslobatto(s)[1]))
    x[begin] = -1; x[end] = 1

    T.(x)
end

lobatto_legendre_points(s; kwargs...) = lobatto_legendre_points(Float64, s; kwargs...)

"""
Lobatto-Legendre nodes on the interval [0,1], i.e., the Lobatto-Legendre points shifted
and scaled from [-1,+1] to [0,1].
"""
function lobatto_legendre_nodes(::Type{T}, s::Integer; IT=BigFloat) where {T}
    T.(shift_nodes(lobatto_legendre_points(IT, s; IT=IT)))
end

lobatto_legendre_nodes(s; kwargs...) = lobatto_legendre_nodes(Float64, s; kwargs...)


function _lobatto_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslobatto(s)
    shift!(b,c)
    QuadratureRule(2s-2, c, b, T)
end


"""
Lobatto-Legendre quadrature.
"""
function LobattoLegendreQuadrature(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    if fast
        return _lobatto_legendre_fast(s, T)
    end

    x = lobatto_legendre_points(IT, s; IT=IT)
    w = [ 2 / ( s*(s-1) * _legendre(s-1, x[i])^2 ) for i in 1:s ]
    return QuadratureRule(2s-2, shift_nodes(x), scale_weights(w), T)
end

LobattoLegendreQuadrature(s; kwargs...) = LobattoLegendreQuadrature(Float64, s; kwargs...)
