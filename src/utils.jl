
# The element types that are computed in BigFloat and rounded at the end, i.e. the numeric
# tower, whose exact members cannot represent the algebraic numbers a quadrature node is in
# general anyway. Everything outside it is taken to compute exactly by itself; see
# _default_arithmetic and _roots.
const _FloatLike = Union{AbstractFloat, Rational, Integer, Complex}

"""
The arithmetic in which nodes and weights are computed by default for the element type `T`.

Types from the numeric tower — floating point, rational, integer and complex ones — are
computed in `BigFloat` and rounded to `T` at the very end, so that the result is accurate to
full precision in `T`. Every other type is computed in itself, on the assumption that it does
its own arithmetic exactly, so that a symbolic `T` yields nodes and weights in closed form.

Should that assumption not hold for some element type, `IT=BigFloat` recovers the behaviour
of the numeric tower for it.
"""
_default_arithmetic(::Type{T}) where {T <: _FloatLike} = BigFloat
_default_arithmetic(::Type{T}) where {T} = T

@doc raw"""
The Legendre polynomial ``P_j(x)`` of degree ``j`` on the interval ``[-1,+1]``, evaluated by
the three-term recurrence

```math
j \, P_j (x) = (2j-1) \, x \, P_{j-1} (x) - (j-1) \, P_{j-2} (x) ,
\qquad P_0 = 1 , \quad P_1 = x .
```

The recurrence is used rather than the equivalent Rodrigues formula

```math
P_j (x) = \frac{1}{j! \, 2^j} \, \frac{d^j}{dx^j} \big( x^2 - 1 \big)^j ,
```

because it costs ``O(j)`` operations — it is carried upwards in a loop — instead of building
and differentiating a polynomial of degree ``2j``. The Rodrigues form is what
identifies the ``(s-2)``-nd derivative of ``(1-x^2)^{s-1}`` used for the Lobatto nodes as an
antiderivative of ``P_{s-1}``; that its roots are the whole Lobatto node set is the separate,
Jacobi-Rodrigues argument given in the manual, cf. [`lobatto_legendre_nodes`](@ref).

Works for any `x` supporting arithmetic, including a `Polynomial` — see
[`QuadratureRules._legendre_polynomial`](@ref) — and symbolic types.
"""
function _legendre(j::Int, x::T) where {T}
    j <= 0 && return one(T)
    j == 1 && return x

    # Carried upwards in a loop rather than by recursing on j-1 and j-2, which would evaluate
    # the same P_k exponentially often. The arithmetic is unchanged, so the result is
    # identical bit for bit.
    p₀, p₁ = one(T), x

    for k in 2:j
        p₀, p₁ = p₁, ((2k-1) * p₁ * x - (k-1) * p₀) / k
    end

    p₁
end

"Legendre polynomial P_s(x) of degree s on the interval [-1..+1]."
function _legendre_polynomial(s, T = BigFloat)
    _legendre(s, Polynomial(T[0, 1]))
end

"""
Refine the approximate roots `x₀` of the polynomial `p` with Newton's method.

The iteration is carried out in the arithmetic of `p`'s coefficients and runs
until the correction stops decreasing, i.e., until the roots are accurate to
the working precision. This is used to compute quadrature nodes in arbitrary
precision from double precision initial guesses. All roots are assumed to be
real and simple.
"""
function _newton_roots(p::Polynomial{T}, x₀::AbstractVector; maxiter = 100) where {T}
    dp = Polynomials.derivative(p)
    x = Vector{T}(undef, length(x₀))

    for i in eachindex(x₀)
        xᵢ = T(x₀[i])
        δ = p(xᵢ) / dp(xᵢ)
        xᵢ -= δ

        for _ in 2:maxiter
            δ̃ = p(xᵢ) / dp(xᵢ)
            xᵢ -= δ̃
            (iszero(δ̃) || abs(δ̃) ≥ abs(δ)) && break
            δ = δ̃
        end

        x[i] = xᵢ
    end

    return x
end

"""
The roots of the polynomial `p`, all assumed to be real and simple.

For coefficients from the numeric tower the double precision approximations `x₀()` are
refined with Newton's method in the arithmetic of `p`. For any other coefficient type, in
particular a symbolic one, the roots are computed exactly as the eigenvalues of the companion
matrix; this needs `eigvals` for that type, which a computer algebra system supplies and which
resolves the roots into radicals. In both cases they come out in no particular order, so the
caller sorts them.

`x₀` is a thunk rather than a vector so that the initial guess, which the exact branch has no
use for, is only computed where it is actually needed.
"""
_roots(p::Polynomial{T}, x₀) where {T <: _FloatLike} = _newton_roots(p, x₀())

_roots(p::Polynomial, x₀) = Polynomials.roots(p)

"Shift and scale nodes from the interval [-1,+1] to the interval [0,1]."
function shift_nodes(c)
    (c .+ 1) ./ 2
end

"Shift and scale nodes from the interval [0,1] to the interval [-1,+1]."
function unshift_nodes(c)
    2 .* c .- 1
end

"Scale weights from the interval [-1,+1] to the interval [0,1]."
function scale_weights(b)
    b ./ 2
end

"Scale weights from the interval [0,1] to the interval [-1,+1]."
function unscale_weights(b)
    b .* 2
end

# Map nodes and weights from the interval their closed form is formulated on to the interval
# the caller asked for. Each family computes on its own interval and converts once at the very
# end, in the working precision, so that the two conventions agree exactly at equal working
# precision and up to rounding across precisions.
_nodes_from_symmetric(x, ::SymmetricInterval) = x
_nodes_from_symmetric(x, ::UnitInterval) = shift_nodes(x)

_nodes_from_unit(c, ::UnitInterval) = c
_nodes_from_unit(c, ::SymmetricInterval) = unshift_nodes(c)

_weights_from_symmetric(w, ::SymmetricInterval) = w
_weights_from_symmetric(w, ::UnitInterval) = scale_weights(w)

_weights_from_unit(b, ::UnitInterval) = b
_weights_from_unit(b, ::SymmetricInterval) = unscale_weights(b)

"Scale nodes and weights from the interval [-1,+1] to the interval [0,1]."
function shift!(b, c)
    b .= b ./ 2
    c .= (c .+ 1) ./ 2
end

"Scale nodes and weights from the interval [0,1] to the interval [-1,+1]."
function unshift!(b, c)
    b .= b .* 2
    c .= 2 .* c .- 1
end
