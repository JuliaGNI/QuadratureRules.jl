# The left Radau-Legendre nodes on [-1,+1], i.e. the roots of P_{s-1} + P_s, one of which is
# set to exactly -1, in their own arithmetic.
function _radau_legendre_nodes(s, ::Val{:left}, IT)
    P = _legendre_polynomial(s-1, IT) + _legendre_polynomial(s, IT)
    x = sort(_roots(P, () -> FastGaussQuadrature.gaussradau(s)[1]))
    x[begin] = -1

    x
end

# The right nodes are the reflection of the left ones. Reflecting rather than evaluating the
# mirrored closed form is what makes the two variants exact mirror images: the recurrence for
# P_{s-1}(-x) does not reproduce P_{s-1}(x) bit for bit.
function _radau_legendre_nodes(s, ::Val{:right}, IT)
    -reverse(_radau_legendre_nodes(s, Val(:left), IT))
end

# The left Radau-Legendre weights on [-1,+1] belonging to the precomputed left nodes `x`, in
# their own arithmetic. Taking the nodes as an argument lets _radau_legendre share the closed
# form with radau_legendre_weights without repeating the root find.
function _radau_legendre_weights(x::AbstractVector{IT}) where {IT}
    s = length(x)

    [(1 - x[i]) / (s^2 * _legendre(s-1, x[i])^2) for i in 1:s]
end

# Nodes and weights on [-1,+1] from a single root find, the right variant obtained by
# reflecting the left one.
function _radau_legendre(s, ::Val{:left}, IT)
    x = _radau_legendre_nodes(s, Val(:left), IT)

    x, _radau_legendre_weights(x)
end

function _radau_legendre(s, ::Val{:right}, IT)
    x, w = _radau_legendre(s, Val(:left), IT)

    -reverse(x), reverse(w)
end

@doc raw"""
    radau_legendre_nodes(s, endpoint; kwargs...)
    radau_legendre_nodes(T, s, ::Val{endpoint}; IT=_default_arithmetic(T),
                         interval=UnitInterval())

The `s` Radau-Legendre nodes, i.e., one prescribed endpoint of the interval together with the
``s-1`` free nodes that maximise the degree of exactness, mapped to `interval`.

Which endpoint is prescribed is selected by `endpoint`:

- `:left` includes the left end of the interval, the classical Gauss-Radau convention and the
  one underlying the Radau IA Runge-Kutta methods,
- `:right` includes the right end, the convention underlying the Radau IIA methods.

There is deliberately no default for `endpoint`: the two variants are not interchangeable,
and silently choosing one would be easy to overlook.

The left nodes are the ``s`` roots of ``P_{s-1} + P_s`` on ``[-1,+1]``, one of which is
exactly ``-1``. They are computed in the arithmetic `IT`, either Newton-refined from the
double precision approximations of `FastGaussQuadrature.gaussradau` or, for a non-numeric
`IT` such as a symbolic one, exactly from the companion matrix, and the prescribed endpoint
is then set to exactly ``-1``. The right nodes are obtained by reflection,
``x \mapsto -x``, which makes the two variants exact mirror images of each other. As the
prescribed endpoint is exact, so is the corresponding node after the mapping.

Note that only the *interval* of `SymmetricInterval` is symmetric: the Radau nodes are
asymmetric by construction on either interval, which is the whole point of the family.

Equivalently, and as the Radau nodes are usually stated in the Runge-Kutta literature, on the
default [`UnitInterval`](@ref) they are the roots of

```math
\frac{d^{\,s-1}}{dx^{\,s-1}} \big( x^s (x - 1)^{s-1} \big)
\qquad \text{for } \texttt{:left} , \qquad
\frac{d^{\,s-1}}{dx^{\,s-1}} \big( x^{s-1} (x - 1)^s \big)
\qquad \text{for } \texttt{:right} ,
```

each of degree ``s``. These agree with ``P_{s-1} + P_s`` through the Rodrigues formula for the
Jacobi polynomial ``P^{(0,1)}_{s-1}``,

```math
P^{(0,1)}_{s-1} (x) \; \propto \; \frac{1}{1+x} \, \frac{d^{\,s-1}}{dx^{\,s-1}}
    \big( (1-x)^{s-1} (1+x)^{s} \big) ,
```

whose *numerator* is what the differentiated product above becomes under the map to
``[-1,+1]``. That numerator is therefore proportional to ``P_{s-1} + P_s`` itself, of degree
``s``, and the division by ``1+x`` is precisely what strips the prescribed endpoint off to
leave the degree ``s-1`` Jacobi polynomial whose roots are the free nodes. The two
expressions are mirror images of one another under ``x \mapsto 1-x``, as the node sets are.

# Arguments
- `T`: element type of the returned vector, `Float64` if omitted.
- `s`: number of nodes.
- `endpoint`: `:left` or `:right`, the endpoint included among the nodes.
- `IT`: arithmetic in which the roots are computed, `BigFloat` for a numeric `T` and `T`
  itself otherwise, cf. [`QuadratureRules._default_arithmetic`](@ref).
- `interval`: [`UnitInterval`](@ref) for ``[0,1]``, the default and the interval of
  [`RadauLegendreQuadrature`](@ref), or [`SymmetricInterval`](@ref) for ``[-1,+1]``.

```jldoctest
julia> radau_legendre_nodes(2, :left)
2-element Vector{Float64}:
 0.0
 0.6666666666666666

julia> radau_legendre_nodes(2, :right)
2-element Vector{Float64}:
 0.3333333333333333
 1.0

julia> radau_legendre_nodes(2, :left; interval = SymmetricInterval())
2-element Vector{Float64}:
 -1.0
  0.3333333333333333
```
"""
function radau_legendre_nodes(::Type{T}, s::Integer, endpoint::Val;
        IT = _default_arithmetic(T),
        interval::QuadratureInterval = UnitInterval()) where {T}
    T.(_nodes_from_symmetric(_radau_legendre_nodes(s, endpoint, IT), interval))
end

function radau_legendre_nodes(s, endpoint; kwargs...)
    radau_legendre_nodes(Float64, s, Val(endpoint); kwargs...)
end

@doc raw"""
    radau_legendre_weights(s, endpoint; kwargs...)
    radau_legendre_weights(T, s, ::Val{endpoint}; IT=_default_arithmetic(T),
                           interval=UnitInterval())

The `s` Radau-Legendre weights belonging to the nodes returned by
[`radau_legendre_nodes`](@ref) for the same `endpoint` and `interval`. All of them are
positive.

They are given in closed form by

```math
w_i = \frac{1 \mp x_i}{s^2 \, \big[ P_{s-1}(x_i) \big]^2} ,
```

with the upper sign for `endpoint = :left` and the lower one for `endpoint = :right`. Like
the corresponding Lobatto formula this holds for the free nodes and for the prescribed
endpoint alike, where ``P_{s-1}(\mp 1)^2 = 1`` reduces it to ``2/s^2``. Being formulated on
``[-1,+1]``, where the weights sum to ``2``, it is the `SymmetricInterval` weights that are
primary here; those on ``[0,1]`` are obtained as ``b_i = w_i / 2`` and sum to ``1``.

In terms of the [`UnitInterval`](@ref) nodes ``c_i`` the same weights read

```math
b_i = \frac{1 \mp (2 c_i - 1)}{2 \, s^2 \, \big[ P_{s-1}(2 c_i - 1) \big]^2} ,
```

with the signs as above, so that the prescribed endpoint again carries ``1/s^2``.

See [`radau_legendre_nodes`](@ref) for the arguments.

```jldoctest
julia> radau_legendre_weights(2, :left)
2-element Vector{Float64}:
 0.25
 0.75

julia> radau_legendre_weights(2, :left; interval = SymmetricInterval())
2-element Vector{Float64}:
 0.5
 1.5
```
"""
function radau_legendre_weights(::Type{T}, s::Integer, endpoint::Val;
        IT = _default_arithmetic(T),
        interval::QuadratureInterval = UnitInterval()) where {T}
    _, w = _radau_legendre(s, endpoint, IT)

    T.(_weights_from_symmetric(w, interval))
end

function radau_legendre_weights(s, endpoint; kwargs...)
    radau_legendre_weights(Float64, s, Val(endpoint); kwargs...)
end

function _radau_legendre_fast(s, T, ::Val{:left})
    c, b = FastGaussQuadrature.gaussradau(s)
    shift!(b, c)
    QuadratureRule(2s-1, c, b, T)
end

function _radau_legendre_fast(s, T, ::Val{:right})
    c, b = FastGaussQuadrature.gaussradau(s)
    c .= -c
    reverse!(c)
    reverse!(b)
    shift!(b, c)
    QuadratureRule(2s-1, c, b, T)
end

@doc raw"""
    RadauLegendreQuadrature(s, endpoint; IT=BigFloat, fast=false)
    RadauLegendreQuadrature(T, s, ::Val{endpoint}; IT=_default_arithmetic(T), fast=false)

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
[`radau_legendre_nodes`](@ref) for `endpoint`.

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
function RadauLegendreQuadrature(::Type{T}, s::Integer, endpoint::Val;
        IT = _default_arithmetic(T), fast = false) where {T}
    if fast
        return _radau_legendre_fast(s, T, endpoint)
    end

    x, w = _radau_legendre(s, endpoint, IT)

    return QuadratureRule(2s-1, shift_nodes(x), scale_weights(w), T)
end

function RadauLegendreQuadrature(s, endpoint; kwargs...)
    RadauLegendreQuadrature(Float64, s, Val(endpoint); kwargs...)
end
