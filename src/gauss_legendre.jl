
_gauss_legendre_points(P, s) = sort(_newton_roots(P, FastGaussQuadrature.gausslegendre(s)[1]))

function _gauss_legendre_points(s::Integer, IT, fast)
    if fast
        return IT.(FastGaussQuadrature.gausslegendre(s)[1])
    end

    _gauss_legendre_points(_legendre_polynomial(s, IT), s)
end

"""
Gauss-Legendre points on the interval [-1,+1].

The points are computed in the arithmetic `IT` and converted to `T`. With `fast = true`
the double precision points of `FastGaussQuadrature` are used instead.
"""
function gauss_legendre_points(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    T.(_gauss_legendre_points(s, IT, fast))
end

gauss_legendre_points(s; kwargs...) = gauss_legendre_points(Float64, s; kwargs...)

"""
Gauss-Legendre nodes on the interval [0,1], i.e., the Gauss-Legendre points shifted
and scaled from [-1,+1] to [0,1].
"""
function gauss_legendre_nodes(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    T.(shift_nodes(_gauss_legendre_points(s, IT, fast)))
end

gauss_legendre_nodes(s; kwargs...) = gauss_legendre_nodes(Float64, s; kwargs...)


function _gauss_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    QuadratureRule(2s, c, b, T)
end

"""
Gauss-Legendre quadrature.
"""
function GaussLegendreQuadrature(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    if fast
        return _gauss_legendre_fast(s, T)
    end

    P = _legendre_polynomial(s, IT)
    D = Polynomials.derivative(P)

    x = _gauss_legendre_points(P, s)

    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(IT[-x[i], 1]) )^2 )
        I(1) - I(-1)
    end

    w = [ inti(i) / D(x[i])^2  for i in 1:s ]

    return QuadratureRule(2s, shift_nodes(x), scale_weights(w), T)
end

GaussLegendreQuadrature(s; kwargs...) = GaussLegendreQuadrature(Float64, s; kwargs...)
