
@doc raw"""
    gauss_legendre_points(s)
    gauss_legendre_points(T, s; IT=BigFloat)

The `s` Gauss-Legendre points on the interval ``[-1,+1]``, i.e., the roots of the
Legendre polynomial ``P_s``.

The roots are obtained by refining the double precision approximations of
`FastGaussQuadrature.gausslegendre` with Newton's method in the arithmetic `IT`, and
are then converted to `T`. Choosing `IT=BigFloat` (the default) therefore yields nodes
that are accurate to full `BigFloat` precision, independently of `T`.

# Arguments
- `T`: element type of the returned vector, `Float64` if omitted.
- `s`: number of points.
- `IT`: arithmetic in which the roots are computed.

```jldoctest
julia> gauss_legendre_points(2)
2-element Vector{Float64}:
 -0.5773502691896257
  0.5773502691896257
```

See also [`gauss_legendre_nodes`](@ref) for the same points on ``[0,1]``.
"""
function gauss_legendre_points(::Type{T}, s::Integer; IT=BigFloat) where {T}
    T.(sort(_newton_roots(_legendre_polynomial(s, IT), FastGaussQuadrature.gausslegendre(s)[1])))
end

gauss_legendre_points(s; kwargs...) = gauss_legendre_points(Float64, s; kwargs...)

@doc raw"""
    gauss_legendre_nodes(s)
    gauss_legendre_nodes(T, s; IT=BigFloat)

The `s` Gauss-Legendre nodes on the interval ``[0,1]``, i.e., the Gauss-Legendre points
shifted and scaled from ``[-1,+1]`` to ``[0,1]`` by ``c_i = (x_i + 1)/2``.

These are the nodes of [`GaussLegendreQuadrature`](@ref). See
[`gauss_legendre_points`](@ref) for the arguments.

```jldoctest
julia> gauss_legendre_nodes(2)
2-element Vector{Float64}:
 0.2113248654051871
 0.7886751345948129
```
"""
function gauss_legendre_nodes(::Type{T}, s::Integer; IT=BigFloat) where {T}
    T.(shift_nodes(gauss_legendre_points(IT, s; IT=IT)))
end

gauss_legendre_nodes(s; kwargs...) = gauss_legendre_nodes(Float64, s; kwargs...)

@doc raw"""
    gauss_legendre_point_weights(s)
    gauss_legendre_point_weights(T, s; IT=BigFloat)

The `s` Gauss-Legendre weights for the interval ``[-1,+1]``, belonging to the points
returned by [`gauss_legendre_points`](@ref) and summing to ``2``.

They are computed from the closed form

```math
w_i = \frac{1}{P_s'(x_i)^2} \int_{-1}^{+1} \left( \frac{P_s(x)}{x - x_i} \right)^2 dx ,
```

where the integrand is the square of the (unnormalised) Lagrange basis polynomial
associated with ``x_i``, evaluated by exact polynomial division and integration in the
arithmetic `IT`. See [`gauss_legendre_points`](@ref) for the arguments.

```jldoctest
julia> gauss_legendre_point_weights(2)
2-element Vector{Float64}:
 1.0
 1.0
```

See also [`gauss_legendre_weights`](@ref) for the same weights on ``[0,1]``.
"""
function gauss_legendre_point_weights(::Type{T}, s::Integer; IT=BigFloat) where {T}
    P = _legendre_polynomial(s, IT)
    D = Polynomials.derivative(P)

    x = gauss_legendre_points(IT, s; IT=IT)

    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(IT[-x[i], 1]) )^2 )
        I(1) - I(-1)
    end

    T.([ inti(i) / D(x[i])^2  for i in 1:s ])
end

gauss_legendre_point_weights(s; kwargs...) = gauss_legendre_point_weights(Float64, s; kwargs...)

@doc raw"""
    gauss_legendre_weights(s)
    gauss_legendre_weights(T, s; IT=BigFloat)

The `s` Gauss-Legendre weights for the interval ``[0,1]``, i.e., the weights of
[`gauss_legendre_point_weights`](@ref) halved so that they sum to ``1``.

These are the weights of [`GaussLegendreQuadrature`](@ref). See
[`gauss_legendre_points`](@ref) for the arguments.

```jldoctest
julia> gauss_legendre_weights(2)
2-element Vector{Float64}:
 0.5
 0.5
```
"""
function gauss_legendre_weights(::Type{T}, s::Integer; IT=BigFloat) where {T}
    T.(scale_weights(gauss_legendre_point_weights(IT, s; IT=IT)))
end

gauss_legendre_weights(s; kwargs...) = gauss_legendre_weights(Float64, s; kwargs...)


function _gauss_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    QuadratureRule(2s, c, b, T)
end

@doc raw"""
    GaussLegendreQuadrature(s; IT=BigFloat, fast=false)
    GaussLegendreQuadrature(T, s; IT=BigFloat, fast=false)

The Gauss-Legendre quadrature rule with `s` nodes on the interval ``[0,1]``.

The nodes are the roots of the Legendre polynomial ``P_s`` mapped to ``[0,1]``. Being
the free choice of both nodes and weights, they achieve the maximal possible degree of
exactness ``2s-1``, so the rule has **order ``2s``**. All nodes lie in the interior of
the interval and all weights are positive.

The weights are computed from the closed form

```math
w_i = \frac{1}{P_s'(x_i)^2} \int_{-1}^{+1} \left( \frac{P_s(x)}{x - x_i} \right)^2 dx ,
```

where the integrand is the square of the (unnormalised) Lagrange basis polynomial
associated with ``x_i``, evaluated by exact polynomial division and integration in the
arithmetic `IT`. The nodes and weights are finally shifted and scaled to ``[0,1]``, and
are also available on their own as [`gauss_legendre_nodes`](@ref) and
[`gauss_legendre_weights`](@ref).

# Arguments
- `T`: element type of the resulting rule, `Float64` if omitted.
- `s`: number of nodes.
- `IT`: arithmetic in which nodes and weights are computed, `BigFloat` by default, so
  that the result is accurate to full precision in `T`.
- `fast`: if `true`, take the nodes and weights directly from
  `FastGaussQuadrature.gausslegendre`, which is much faster but computes in double
  precision, so the result is only accurate to about `Float64` precision.

```jldoctest
julia> GaussLegendreQuadrature(2)
QuadratureRule{Float64, 2}(4, [0.2113248654051871, 0.7886751345948129], [0.5, 0.5])

julia> GaussLegendreQuadrature(3)(x -> x^5)   # exact up to degree 5
0.16666666666666669
```

See also [`LobattoLegendreQuadrature`](@ref), [`gauss_legendre_nodes`](@ref).
"""
function GaussLegendreQuadrature(::Type{T}, s::Integer; IT=BigFloat, fast=false) where {T}
    if fast
        return _gauss_legendre_fast(s, T)
    end

    c = gauss_legendre_nodes(IT, s; IT=IT)
    b = gauss_legendre_weights(IT, s; IT=IT)

    return QuadratureRule(2s, c, b, T)
end

GaussLegendreQuadrature(s; kwargs...) = GaussLegendreQuadrature(Float64, s; kwargs...)
