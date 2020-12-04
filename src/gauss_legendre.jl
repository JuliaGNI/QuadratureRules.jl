
function _gauss_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    QuadratureRule(2s, c, b, T)
end

function _ShiftedLegendre(s, T=BigFloat)
    p = vcat([zero(T) for _ in 1:s], one(T))
    ShiftedLegendre(p)
end

function GaussLegendreQuadrature(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    if fast
        return _gauss_legendre_fast(s, T)
    end

    P = convert(Polynomial, _ShiftedLegendre(s,IT))
    D = Polynomials.derivative(P)

    c = sort(IT.(Polynomials.roots(P)))

    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(IT[-c[i], 1]) )^2 )
        I(1) - I(0)
    end
    
    b = [ inti(i) / D(c[i])^2  for i in 1:s ]

    return QuadratureRule(2s, c, b, T)
end

GaussLegendreQuadrature(s; kwargs...) = GaussLegendreQuadrature(Float64, s; kwargs...)
