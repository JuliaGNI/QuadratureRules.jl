# The Lobatto-Legendre nodes on [-1,+1], i.e. the endpoints together with the roots of
# P'_{s-1}, in their own arithmetic. Rather than differentiating the Legendre polynomial, the
# interior nodes are obtained as the roots of the (s-2)-nd derivative of (1-x^2)^{s-1}, which
# has the same roots by Rodrigues' formula.
function _lobatto_legendre_nodes(s, IT)
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    D = Polynomials.derivative(Polynomial(IT[1, 0, -1])^(s-1), s-2)
    x = sort(_roots(D, () -> FastGaussQuadrature.gausslobatto(s)[1]))
    x[begin] = -1; x[end] = 1

    x
end

@doc raw"""
    lobatto_legendre_nodes(s; kwargs...)
    lobatto_legendre_nodes(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Lobatto-Legendre nodes, i.e., the two endpoints of the interval together with the
``s-2`` roots of ``P_{s-1}'`` mapped to `interval`.

Rather than differentiating the Legendre polynomial, the interior nodes are obtained as the
roots of the ``(s-2)``-nd derivative of ``(1-x^2)^{s-1}``, which has the same roots by
Rodrigues' formula. That derivative has degree ``s``, and its ``s`` roots are therefore the
whole node set, endpoints included. The roots are computed in the arithmetic `IT`, either
Newton-refined from the double precision approximations of
`FastGaussQuadrature.gausslobatto` or, for a non-numeric `IT` such as a symbolic one, exactly
from the companion matrix.

On the default [`UnitInterval`](@ref) the corresponding polynomial is

```math
\frac{d^{\,s-2}}{dx^{\,s-2}} \big( (x - x^2)^{s-1} \big) ,
```

since ``x - x^2 = (1 - \xi^2)/4`` under ``x = (\xi+1)/2``; this is the form in which the
Lobatto nodes are usually stated in the Runge-Kutta literature.

The endpoints are set to exactly ``\mp 1`` before the mapping, so the first and last node
come out as exactly `-1` and `+1` on the symmetric interval and exactly `0` and `1` on the
unit one.

Throws an `ErrorException` for `s == 1`, as a Lobatto rule needs at least the two
endpoints. See [`gauss_legendre_nodes`](@ref) for the arguments.

```jldoctest
julia> lobatto_legendre_nodes(3)
3-element Vector{Float64}:
 0.0
 0.5
 1.0

julia> lobatto_legendre_nodes(3; interval = SymmetricInterval())
3-element Vector{Float64}:
 -1.0
  0.0
  1.0
```
"""
function lobatto_legendre_nodes(::Type{T}, s::Integer; IT=_default_arithmetic(T),
                                interval::QuadratureInterval=UnitInterval()) where {T}
    T.(_nodes_from_symmetric(_lobatto_legendre_nodes(s, IT), interval))
end

lobatto_legendre_nodes(s; kwargs...) = lobatto_legendre_nodes(Float64, s; kwargs...)

# The Lobatto-Legendre weights on [-1,+1] belonging to the precomputed nodes `x`, in their
# own arithmetic. Taking the nodes as an argument lets the constructor share the closed
# form with lobatto_legendre_weights without repeating the root find.
function _lobatto_legendre_weights(x::AbstractVector{IT}) where {IT}
    s = length(x)

    [ 2 / ( s*(s-1) * _legendre(s-1, x[i])^2 )  for i in 1:s ]
end

@doc raw"""
    lobatto_legendre_weights(s; kwargs...)
    lobatto_legendre_weights(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Lobatto-Legendre weights belonging to the nodes returned by
[`lobatto_legendre_nodes`](@ref) for the same `interval`. All of them are positive.

They are given in closed form by

```math
w_i = \frac{2}{s \, (s-1) \, \big[ P_{s-1}(x_i) \big]^2} ,
```

which holds for the interior nodes and for the two endpoints alike. Being formulated on
``[-1,+1]``, where the weights sum to ``2``, it is the `SymmetricInterval` weights that are
primary here; those on ``[0,1]`` are obtained as ``b_i = w_i / 2`` and sum to ``1``.

In terms of the [`UnitInterval`](@ref) nodes ``c_j`` the halving cancels the numerator, so the
same weights read

```math
b_j = \frac{1}{s \, (s-1) \, \big[ P_{s-1}(2 c_j - 1) \big]^2} ,
```

which is the form in which they are usually tabulated for the Lobatto Runge-Kutta methods.

Throws an `ErrorException` for `s == 1`. See [`gauss_legendre_nodes`](@ref) for the
arguments.

```jldoctest
julia> lobatto_legendre_weights(3)
3-element Vector{Float64}:
 0.16666666666666666
 0.6666666666666666
 0.16666666666666666

julia> lobatto_legendre_weights(3; interval = SymmetricInterval())
3-element Vector{Float64}:
 0.3333333333333333
 1.3333333333333333
 0.3333333333333333
```
"""
function lobatto_legendre_weights(::Type{T}, s::Integer; IT=_default_arithmetic(T),
                                  interval::QuadratureInterval=UnitInterval()) where {T}
    w = _lobatto_legendre_weights(_lobatto_legendre_nodes(s, IT))

    T.(_weights_from_symmetric(w, interval))
end

lobatto_legendre_weights(s; kwargs...) = lobatto_legendre_weights(Float64, s; kwargs...)


function _lobatto_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslobatto(s)
    shift!(b,c)
    QuadratureRule(2s-2, c, b, T)
end


@doc raw"""
    LobattoLegendreQuadrature(s; IT=BigFloat, fast=false)
    LobattoLegendreQuadrature(T, s; IT=_default_arithmetic(T), fast=false)

The Lobatto-Legendre quadrature rule with `s` nodes on the interval ``[0,1]``.

In contrast to [`GaussLegendreQuadrature`](@ref), both endpoints of the interval are
included among the nodes. This constraint costs two degrees of exactness: only the
``s-2`` interior nodes are free, so the maximal degree of exactness is ``2s-3`` and the
rule has **order ``2s-2``**. In exchange, the rule can be applied when the values at the
endpoints are needed anyway, which is why it is common in finite element methods,
collocation schemes and variational integrators.

The interior nodes are the roots of ``P_{s-1}'`` and the weights are given in closed
form by

```math
w_i = \frac{2}{s \, (s-1) \, \big[ P_{s-1}(x_i) \big]^2} ,
```

which holds for the endpoints as well. Nodes and weights are computed in the arithmetic
`IT` and then shifted and scaled to ``[0,1]``, and are also available on their own as
[`lobatto_legendre_nodes`](@ref) and [`lobatto_legendre_weights`](@ref).

Throws an `ErrorException` for `s == 1`. See [`GaussLegendreQuadrature`](@ref) for the
meaning of `T`, `IT` and `fast`.

```jldoctest
julia> LobattoLegendreQuadrature(2) == TrapezoidalQuadrature()
true

julia> LobattoLegendreQuadrature(3)
QuadratureRule{Float64, 3}(4, [0.0, 0.5, 1.0], [0.16666666666666666, 0.6666666666666666, 0.16666666666666666])
```

The three-node rule is Simpson's rule, and it is exact up to degree ``2 \cdot 3 - 3 = 3``.
"""
function LobattoLegendreQuadrature(::Type{T}, s::Integer; IT=_default_arithmetic(T), fast=false) where {T}
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    if fast
        return _lobatto_legendre_fast(s, T)
    end

    x = _lobatto_legendre_nodes(s, IT)
    w = _lobatto_legendre_weights(x)

    return QuadratureRule(2s-2, shift_nodes(x), scale_weights(w), T)
end

LobattoLegendreQuadrature(s; kwargs...) = LobattoLegendreQuadrature(Float64, s; kwargs...)
