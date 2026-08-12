
@doc raw"""
    radau_legendre_points(s, endpoint; IT)
    radau_legendre_points(T, s, ::Val{endpoint}; IT)

The `s` Radau-Legendre points on the interval ``[-1,+1]``, i.e., one prescribed endpoint
together with the ``s-1`` free points that maximise the degree of exactness.

Which endpoint is prescribed is selected by `endpoint`:

- `:left` includes ``-1``, the classical Gauss-Radau convention and the one underlying the
  Radau IA Runge-Kutta methods,
- `:right` includes ``+1``, the convention underlying the Radau IIA methods.

The left points are the ``s`` roots of ``P_{s-1} + P_s``, one of which is exactly ``-1``.
They are computed in the arithmetic `IT`, either Newton-refined from the double precision
approximations of `FastGaussQuadrature.gaussradau` or, for an exact `IT` such as a symbolic
one, exactly from the companion matrix, and the prescribed endpoint is then set to exactly
``-1``. The right points are obtained by reflection, ``x \mapsto -x``, which makes the two
variants exact mirror images of each other.

There is deliberately no default for `endpoint`: the two variants are not interchangeable,
and silently choosing one would be easy to overlook.

# Arguments
- `T`: element type of the returned vector, `Float64` if omitted.
- `s`: number of points.
- `endpoint`: `:left` or `:right`, the endpoint included among the points.
- `IT`: arithmetic in which the roots are computed, `BigFloat` for a floating point `T` and
  `T` itself otherwise.

```jldoctest
julia> radau_legendre_points(2, :left)
2-element Vector{Float64}:
 -1.0
  0.3333333333333333

julia> radau_legendre_points(2, :right)
2-element Vector{Float64}:
 -0.3333333333333333
  1.0
```

See also [`radau_legendre_nodes`](@ref) for the same points on ``[0,1]``.
"""
function radau_legendre_points(::Type{T}, s::Integer, ::Val{:left}; IT=_default_arithmetic(T)) where {T}
    P = _legendre_polynomial(s-1, IT) + _legendre_polynomial(s, IT)
    x = sort(_roots(P, FastGaussQuadrature.gaussradau(s)[1]))
    x[begin] = -1

    T.(x)
end

function radau_legendre_points(::Type{T}, s::Integer, ::Val{:right}; IT=_default_arithmetic(T)) where {T}
    T.(-reverse(radau_legendre_points(IT, s, Val(:left); IT=IT)))
end

radau_legendre_points(s, endpoint; kwargs...) = radau_legendre_points(Float64, s, Val(endpoint); kwargs...)

@doc raw"""
    radau_legendre_nodes(s, endpoint; IT)
    radau_legendre_nodes(T, s, ::Val{endpoint}; IT)

The `s` Radau-Legendre nodes on the interval ``[0,1]``, i.e., the Radau-Legendre points
shifted and scaled from ``[-1,+1]`` to ``[0,1]``.

As the prescribed endpoint of [`radau_legendre_points`](@ref) is exact, the first node is
exactly `0` for `endpoint = :left` and the last node exactly `1` for `endpoint = :right`.
These are the nodes of [`RadauLegendreQuadrature`](@ref). See
[`radau_legendre_points`](@ref) for the arguments.

```jldoctest
julia> radau_legendre_nodes(2, :left)
2-element Vector{Float64}:
 0.0
 0.6666666666666666

julia> radau_legendre_nodes(2, :right)
2-element Vector{Float64}:
 0.3333333333333333
 1.0
```
"""
function radau_legendre_nodes(::Type{T}, s::Integer, endpoint::Val; IT=_default_arithmetic(T)) where {T}
    T.(shift_nodes(radau_legendre_points(IT, s, endpoint; IT=IT)))
end

radau_legendre_nodes(s, endpoint; kwargs...) = radau_legendre_nodes(Float64, s, Val(endpoint); kwargs...)

# The left Radau-Legendre weights on [-1,+1] belonging to the precomputed left points `x`,
# in their own arithmetic. Taking the points as an argument lets _radau_legendre share the
# closed form with radau_legendre_point_weights without repeating the root find.
function _radau_legendre_point_weights(x::AbstractVector{IT}) where {IT}
    s = length(x)

    [ (1 - x[i]) / ( s^2 * _legendre(s-1, x[i])^2 )  for i in 1:s ]
end

# Points and weights on [-1,+1] from a single root find, the right variant obtained by
# reflecting the left one. Reflecting rather than evaluating the mirrored closed form is
# what makes the two variants exact mirror images: the recurrence for P_{s-1}(-x) does not
# reproduce P_{s-1}(x) bit for bit.
function _radau_legendre(s, ::Val{:left}, IT)
    x = radau_legendre_points(IT, s, Val(:left); IT=IT)

    x, _radau_legendre_point_weights(x)
end

function _radau_legendre(s, ::Val{:right}, IT)
    x, w = _radau_legendre(s, Val(:left), IT)

    -reverse(x), reverse(w)
end

@doc raw"""
    radau_legendre_point_weights(s, endpoint; IT)
    radau_legendre_point_weights(T, s, ::Val{endpoint}; IT)

The `s` Radau-Legendre weights for the interval ``[-1,+1]``, belonging to the points
returned by [`radau_legendre_points`](@ref) and summing to ``2``.

They are given in closed form by

```math
w_i = \frac{1 \mp x_i}{s^2 \, \big[ P_{s-1}(x_i) \big]^2} ,
```

with the upper sign for `endpoint = :left` and the lower one for `endpoint = :right`. Like
the corresponding Lobatto formula this holds for the free points and for the prescribed
endpoint alike, where ``P_{s-1}(\mp 1)^2 = 1`` reduces it to ``2/s^2``.

See [`radau_legendre_points`](@ref) for the arguments.

```jldoctest
julia> radau_legendre_point_weights(2, :left)
2-element Vector{Float64}:
 0.5
 1.5
```

See also [`radau_legendre_weights`](@ref) for the same weights on ``[0,1]``.
"""
function radau_legendre_point_weights(::Type{T}, s::Integer, endpoint::Val; IT=_default_arithmetic(T)) where {T}
    _, w = _radau_legendre(s, endpoint, IT)

    T.(w)
end

radau_legendre_point_weights(s, endpoint; kwargs...) = radau_legendre_point_weights(Float64, s, Val(endpoint); kwargs...)

@doc raw"""
    radau_legendre_weights(s, endpoint; IT)
    radau_legendre_weights(T, s, ::Val{endpoint}; IT)

The `s` Radau-Legendre weights for the interval ``[0,1]``, i.e., the weights of
[`radau_legendre_point_weights`](@ref) halved so that they sum to ``1``.

These are the weights of [`RadauLegendreQuadrature`](@ref). See
[`radau_legendre_points`](@ref) for the arguments.

```jldoctest
julia> radau_legendre_weights(2, :left)
2-element Vector{Float64}:
 0.25
 0.75
```
"""
function radau_legendre_weights(::Type{T}, s::Integer, endpoint::Val; IT=_default_arithmetic(T)) where {T}
    T.(scale_weights(radau_legendre_point_weights(IT, s, endpoint; IT=IT)))
end

radau_legendre_weights(s, endpoint; kwargs...) = radau_legendre_weights(Float64, s, Val(endpoint); kwargs...)


function _radau_legendre_fast(s, T, ::Val{:left})
    c, b = FastGaussQuadrature.gaussradau(s)
    shift!(b,c)
    QuadratureRule(2s-1, c, b, T)
end

function _radau_legendre_fast(s, T, ::Val{:right})
    c, b = FastGaussQuadrature.gaussradau(s)
    c .= -c
    reverse!(c); reverse!(b)
    shift!(b,c)
    QuadratureRule(2s-1, c, b, T)
end


@doc raw"""
    RadauLegendreQuadrature(s, endpoint; IT, fast=false)
    RadauLegendreQuadrature(T, s, ::Val{endpoint}; IT, fast=false)

The Radau-Legendre quadrature rule with `s` nodes on the interval ``[0,1]``.

It sits between [`GaussLegendreQuadrature`](@ref), which prescribes no node, and
[`LobattoLegendreQuadrature`](@ref), which prescribes both endpoints: exactly one endpoint
is included among the nodes, selected by `endpoint` as `:left` (the node `0`) or `:right`
(the node `1`). This single constraint costs one degree of exactness, so the maximal degree
of exactness is ``2s-2`` and the rule has **order ``2s-1``**. All weights are positive.

The one-sided constraint is what makes the Radau nodes the natural choice for stiffly
accurate collocation: the Radau IA and Radau IIA Runge-Kutta methods are built on the
`:left` and `:right` nodes respectively.

The free nodes are the roots of ``P_{s-1} + P_s`` other than the prescribed endpoint, and
the weights are given in closed form by

```math
w_i = \frac{1 \mp x_i}{s^2 \, \big[ P_{s-1}(x_i) \big]^2} ,
```

with the upper sign for `:left` and the lower one for `:right`, valid at the prescribed
endpoint as well. Nodes and weights are computed in the arithmetic `IT` and then shifted
and scaled to ``[0,1]``, and are also available on their own as
[`radau_legendre_nodes`](@ref) and [`radau_legendre_weights`](@ref).

Unlike the Lobatto rules, the Radau rules are defined for `s == 1`, where the single node
is the prescribed endpoint and the rule reduces to a Riemann sum.

See [`GaussLegendreQuadrature`](@ref) for the meaning of `T`, `IT` and `fast`, and
[`radau_legendre_points`](@ref) for `endpoint`.

```jldoctest
julia> RadauLegendreQuadrature(1, :left) == RiemannQuadratureLeft()
true

julia> RadauLegendreQuadrature(1, :right) == RiemannQuadratureRight()
true

julia> RadauLegendreQuadrature(2, :right)
QuadratureRule{Float64, 2}(3, [0.3333333333333333, 1.0], [0.75, 0.25])

julia> RadauLegendreQuadrature(3, :right)(x -> x^4)   # exact up to degree 2*3-2
0.19999999999999996
```

See also [`GaussLegendreQuadrature`](@ref) and [`LobattoLegendreQuadrature`](@ref).
"""
function RadauLegendreQuadrature(::Type{T}, s::Integer, endpoint::Val; IT=_default_arithmetic(T), fast=false) where {T}
    if fast
        return _radau_legendre_fast(s, T, endpoint)
    end

    x, w = _radau_legendre(s, endpoint, IT)

    return QuadratureRule(2s-1, shift_nodes(x), scale_weights(w), T)
end

RadauLegendreQuadrature(s, endpoint; kwargs...) = RadauLegendreQuadrature(Float64, s, Val(endpoint); kwargs...)
