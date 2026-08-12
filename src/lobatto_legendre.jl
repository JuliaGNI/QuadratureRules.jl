
@doc raw"""
    lobatto_legendre_points(s)
    lobatto_legendre_points(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Legendre points on the interval ``[-1,+1]``, i.e., the two endpoints
``\pm 1`` together with the ``s-2`` roots of ``P_{s-1}'``.

Rather than differentiating the Legendre polynomial, the interior points are obtained
as the roots of the ``(s-2)``-nd derivative of ``(1-x^2)^{s-1}``, which has the same
roots by Rodrigues' formula. The roots are computed in the arithmetic `IT`, either
Newton-refined from the double precision approximations of
`FastGaussQuadrature.gausslobatto` or, for a non-numeric `IT` such as a symbolic one,
exactly from the companion matrix, and the endpoints are set to exactly ``\pm 1``.

Throws an `ErrorException` for `s == 1`, as a Lobatto rule needs at least the two
endpoints. See [`gauss_legendre_points`](@ref) for the arguments.

```jldoctest
julia> lobatto_legendre_points(3)
3-element Vector{Float64}:
 -1.0
  0.0
  1.0
```
"""
function lobatto_legendre_points(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    if s == 1
        throw(ErrorException("Lobatto quadrature is not defined for one stage."))
    end

    D = Polynomials.derivative(Polynomial(IT[1, 0, -1])^(s-1), s-2)
    x = sort(_roots(D, () -> FastGaussQuadrature.gausslobatto(s)[1]))
    x[begin] = -1; x[end] = 1

    T.(x)
end

lobatto_legendre_points(s; kwargs...) = lobatto_legendre_points(Float64, s; kwargs...)

@doc raw"""
    lobatto_legendre_nodes(s)
    lobatto_legendre_nodes(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Legendre nodes on the interval ``[0,1]``, i.e., the Lobatto-Legendre
points shifted and scaled from ``[-1,+1]`` to ``[0,1]``.

As the endpoints of [`lobatto_legendre_points`](@ref) are exact, the first and last
node are exactly `0` and `1`. These are the nodes of
[`LobattoLegendreQuadrature`](@ref).

```jldoctest
julia> lobatto_legendre_nodes(3)
3-element Vector{Float64}:
 0.0
 0.5
 1.0
```
"""
function lobatto_legendre_nodes(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    T.(shift_nodes(lobatto_legendre_points(IT, s; IT=IT)))
end

lobatto_legendre_nodes(s; kwargs...) = lobatto_legendre_nodes(Float64, s; kwargs...)

# The Lobatto-Legendre weights on [-1,+1] belonging to the precomputed points `x`, in their
# own arithmetic. Taking the points as an argument lets the constructor share the closed
# form with lobatto_legendre_point_weights without repeating the root find.
function _lobatto_legendre_point_weights(x::AbstractVector{IT}) where {IT}
    s = length(x)

    [ 2 / ( s*(s-1) * _legendre(s-1, x[i])^2 )  for i in 1:s ]
end

@doc raw"""
    lobatto_legendre_point_weights(s)
    lobatto_legendre_point_weights(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Legendre weights for the interval ``[-1,+1]``, belonging to the points
returned by [`lobatto_legendre_points`](@ref) and summing to ``2``.

They are given in closed form by

```math
w_i = \frac{2}{s \, (s-1) \, \big[ P_{s-1}(x_i) \big]^2} ,
```

which holds for the interior points and for the two endpoints alike. Throws an
`ErrorException` for `s == 1`. See [`lobatto_legendre_points`](@ref) for the arguments.

```jldoctest
julia> lobatto_legendre_point_weights(3)
3-element Vector{Float64}:
 0.3333333333333333
 1.3333333333333333
 0.3333333333333333
```

See also [`lobatto_legendre_weights`](@ref) for the same weights on ``[0,1]``.
"""
function lobatto_legendre_point_weights(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    T.(_lobatto_legendre_point_weights(lobatto_legendre_points(IT, s; IT=IT)))
end

lobatto_legendre_point_weights(s; kwargs...) = lobatto_legendre_point_weights(Float64, s; kwargs...)

@doc raw"""
    lobatto_legendre_weights(s)
    lobatto_legendre_weights(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Legendre weights for the interval ``[0,1]``, i.e., the weights of
[`lobatto_legendre_point_weights`](@ref) halved so that they sum to ``1``.

These are the weights of [`LobattoLegendreQuadrature`](@ref). Throws an `ErrorException`
for `s == 1`. See [`lobatto_legendre_points`](@ref) for the arguments.

```jldoctest
julia> lobatto_legendre_weights(3)
3-element Vector{Float64}:
 0.16666666666666666
 0.6666666666666666
 0.16666666666666666
```
"""
function lobatto_legendre_weights(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    T.(scale_weights(lobatto_legendre_point_weights(IT, s; IT=IT)))
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

    x = lobatto_legendre_points(IT, s; IT=IT)
    w = _lobatto_legendre_point_weights(x)

    return QuadratureRule(2s-2, shift_nodes(x), scale_weights(w), T)
end

LobattoLegendreQuadrature(s; kwargs...) = LobattoLegendreQuadrature(Float64, s; kwargs...)
