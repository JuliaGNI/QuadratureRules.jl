"""
    clenshaw_curtis_points(s; IT=BigFloat)
    clenshaw_curtis_points(T, s; IT=BigFloat)

The `s` Clenshaw-Curtis points on the interval ``[-1,+1]``, i.e., the Chebyshev points
of the second kind, cf. [`chebyshev_points`](@ref). They include both endpoints.

Identical to [`lobatto_chebyshev_points`](@ref). Requires `s ≥ 2`.
"""
clenshaw_curtis_points(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_points(T, s, Val(2); kwargs...)
clenshaw_curtis_points(s; kwargs...) = clenshaw_curtis_points(Float64, s; kwargs...)

"""
    clenshaw_curtis_nodes(s; IT=BigFloat)
    clenshaw_curtis_nodes(T, s; IT=BigFloat)

The `s` Clenshaw-Curtis nodes on the interval ``[0,1]``, i.e., the Chebyshev nodes of
the second kind, cf. [`chebyshev_nodes`](@ref). The first and last node are `0` and `1`.

These are exactly the nodes of [`ClenshawCurtisQuadrature`](@ref) at the same `T` and `IT`.
Requires `s ≥ 2`.

```jldoctest
julia> clenshaw_curtis_nodes(5)
5-element Vector{Float64}:
 0.0
 0.14644660940672624
 0.5
 0.8535533905932737
 1.0
```
"""
clenshaw_curtis_nodes(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_nodes(T, s, Val(2); kwargs...)
clenshaw_curtis_nodes(s; kwargs...) = clenshaw_curtis_nodes(Float64, s; kwargs...)

@doc raw"""
    clenshaw_curtis_weights(s; IT=BigFloat)
    clenshaw_curtis_weights(T, s; IT=BigFloat)

The `s` Clenshaw-Curtis weights for the interval ``[0,1]``, belonging to the nodes
returned by [`clenshaw_curtis_nodes`](@ref) and summing to ``1``.

They are the weights of the interpolatory rule through the Chebyshev points of the second
kind, given explicitly by

```math
w_k = \frac{c_k}{2n} \left( 1 - \sum_{j=1}^{\lfloor n/2 \rfloor}
      \frac{b_j}{4 j^2 - 1} \cos ( j \vartheta_k ) \right) ,
      \qquad \vartheta_k = \frac{2 \pi k}{n} , \qquad n = s-1 ,
```

with ``c_k = 1`` at the endpoints and ``2`` otherwise, and ``b_j = 1`` for the last term
of an even sum and ``2`` otherwise. All weights are positive. See
[`ClenshawCurtisQuadrature`](@ref) for the derivation, the literature and the arguments.

Throws an `ErrorException` for `s == 1`.

```jldoctest
julia> clenshaw_curtis_weights(3)
3-element Vector{Float64}:
 0.16666666666666666
 0.6666666666666666
 0.16666666666666666
```

See also [`clenshaw_curtis_point_weights`](@ref) for the same weights on ``[-1,+1]``.
"""
function clenshaw_curtis_weights(::Type{T}, s::Integer; IT=BigFloat) where {T}
    if s == 1
        throw(ErrorException("Clenshaw-Curtis quadrature is not defined for one stage."))
    end

    b(j,n) = fld(n,2) == j == cld(n,2) ? IT(1) : IT(2)
    c(k,n) = mod(k,n) == 0 ? IT(1) : IT(2)
    ϑ(k,n) = 2 * IT(π) * k / n

    cctrm(k,j,n) = b(j,n) / IT( 4 * j^2 - 1 ) * cos(j * ϑ(k,n))
    ccsum(k,n) = c(k,n) / IT(2n) * ( IT(1) - mapreduce(j -> cctrm(k,j,n), +, 1:div(n,2); init = zero(IT)) )

    T.([ccsum(i-1,s-1) for i in 1:s])
end

clenshaw_curtis_weights(s; kwargs...) = clenshaw_curtis_weights(Float64, s; kwargs...)

@doc raw"""
    clenshaw_curtis_point_weights(s; IT=BigFloat)
    clenshaw_curtis_point_weights(T, s; IT=BigFloat)

The `s` Clenshaw-Curtis weights for the interval ``[-1,+1]``, i.e., the weights of
[`clenshaw_curtis_weights`](@ref) doubled so that they sum to ``2``.

Throws an `ErrorException` for `s == 1`.

```jldoctest
julia> clenshaw_curtis_point_weights(3)
3-element Vector{Float64}:
 0.3333333333333333
 1.3333333333333333
 0.3333333333333333
```
"""
function clenshaw_curtis_point_weights(::Type{T}, s::Integer; IT=BigFloat) where {T}
    T.(unscale_weights(clenshaw_curtis_weights(IT, s; IT=IT)))
end

clenshaw_curtis_point_weights(s; kwargs...) = clenshaw_curtis_point_weights(Float64, s; kwargs...)


@doc raw"""
    ClenshawCurtisQuadrature(s; IT=BigFloat)
    ClenshawCurtisQuadrature(T, s; IT=BigFloat)

The Clenshaw-Curtis quadrature rule with `s` nodes on the interval ``[0,1]``.

The nodes are the Chebyshev points of the second kind, which include both endpoints of
the interval. The weights are those of the interpolatory rule through them: the
integrand is expanded in a Chebyshev series, whose coefficients follow from the values
at the nodes by a discrete cosine transform, and the series is integrated term by term.
Since ``\int_{-1}^{+1} T_{2j}(x) \, dx = -2/(4j^2-1)`` and the odd terms integrate to
zero, this yields

```math
w_k = \frac{c_k}{n} \left( 1 - \sum_{j=1}^{\lfloor n/2 \rfloor}
      \frac{b_j}{4 j^2 - 1} \cos ( j \vartheta_k ) \right) ,
      \qquad \vartheta_k = \frac{2 \pi k}{n} , \qquad n = s-1 ,
```

with ``c_k = 1`` at the endpoints and ``2`` otherwise, and ``b_j = 1`` for the last term
of an even sum and ``2`` otherwise; a further factor ``1/2`` maps the rule to ``[0,1]``.
This is the explicit form given by Reid; see the
[Clenshaw-Curtis](@ref "Clenshaw-Curtis quadrature") section of the manual for the
derivation and for the literature. Nodes and weights are also available on their own as
[`clenshaw_curtis_nodes`](@ref) and [`clenshaw_curtis_weights`](@ref).

Being an interpolatory rule on `s` nodes, it is exact for polynomials of degree
``\le s-1``. For odd `s` it gains one further degree, because the monomial of degree `s`
that it would otherwise fail on is odd about the midpoint of the interval and its error
therefore cancels. The **order is `s` for even `s` and `s+1` for odd `s`**, and it is
sharp: the rule is not exact one degree beyond. All weights are positive, which is what
guarantees convergence for every continuous integrand.

Although its order is only about half that of [`GaussLegendreQuadrature`](@ref) with the
same number of nodes, Clenshaw-Curtis converges at essentially the same rate for
integrands that are not analytic in a sizable neighbourhood of the interval, and its
nodes are given in closed form, requiring no root finding.

Throws an `ErrorException` for `s == 1`.

# Arguments
- `T`: element type of the resulting rule, `Float64` if omitted.
- `s`: number of nodes.
- `IT`: arithmetic in which nodes and weights are computed, `BigFloat` by default.

!!! note "Why `BigFloat` is the default"
    The weight sum costs ``O(s^2)`` operations, so a lower working precision such as
    `IT=Float64` is considerably faster. It is not the default, however, because the sum
    accumulates round-off in its intermediate terms: computing in `BigFloat` and rounding
    only the final result guarantees weights that are correct to the full precision of `T`,
    whereas `IT=Float64` merely gets close to it. Lower the working precision only when the
    cost matters and a few units in the last place do not.

```jldoctest
julia> ClenshawCurtisQuadrature(3)     # Simpson's rule
QuadratureRule{Float64, 3}(4, [0.0, 0.5, 1.0], [0.16666666666666666, 0.6666666666666666, 0.16666666666666666])

julia> ClenshawCurtisQuadrature(3) == LobattoLegendreQuadrature(3)
true

julia> ClenshawCurtisQuadrature(5)(x -> x^4)
0.19999999999999996
```

See also [`LobattoChebyshevQuadrature`](@ref), which is the same rule, and
[`GaussChebyshevQuadrature`](@ref), which is the analogous rule on the Chebyshev points
of the first kind.
"""
function ClenshawCurtisQuadrature(::Type{T}, s::Integer; IT=BigFloat) where {T}
    if s == 1
        throw(ErrorException("Clenshaw-Curtis quadrature is not defined for one stage."))
    end

    x = clenshaw_curtis_nodes(IT, s; IT=IT)
    w = clenshaw_curtis_weights(IT, s; IT=IT)

    QuadratureRule(isodd(s) ? s+1 : s, x, w, T)
end

ClenshawCurtisQuadrature(s; kwargs...) = ClenshawCurtisQuadrature(Float64, s; kwargs...)
