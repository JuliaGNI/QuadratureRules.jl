
function _lobatto_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslobatto(s)
    shift!(b,c)
    QuadratureRule(2s-2, c, b, T)
end


function LobattoLegendreQuadrature(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    if fast
        return _lobatto_legendre_fast(s, T)
    end

    D(k) = Polynomials.derivative(Polynomial(IT[0, 1, -1])^(k-1), k-2)
    P(k,x) = Polynomials.derivative(Polynomial(IT[-1, 0, 1])^k, k)(x) / factorial(k) / 2^k
    c = sort(IT.(Polynomials.roots(D(s)))); c[begin] = 0; c[end] = 1;
    b = [ 1 / ( s*(s-1) * P(s-1, 2c[i] - 1)^2 ) for i in 1:s ]
    return QuadratureRule(2s-2, c, b, T)
end

LobattoLegendreQuadrature(s; kwargs...) = LobattoLegendreQuadrature(Float64, s; kwargs...)
