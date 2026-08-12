
@doc raw"""
    chebyshev_points(s, kind; IT=BigFloat)
    chebyshev_points(T, s, ::Val{kind}; IT=_default_arithmetic(T))

The `s` Chebyshev points of the first (`kind = 1`) or second (`kind = 2`) kind on the
interval ``[-1,+1]``, in ascending order.

The points of the first kind are the roots of the Chebyshev polynomial ``T_s``,

```math
x_i = \cos \left( \frac{(2i-1) \, \pi}{2s} \right) , \qquad i = 1, \dots, s ,
```

which lie strictly inside the interval. They are evaluated in the equivalent `sin` form
used in the implementation, which is more accurate near the ends of the interval. The
points of the second kind are the extrema of ``T_{s-1}``,

```math
x_i = \cos \left( \frac{(i-1) \, \pi}{s-1} \right) , \qquad i = 1, \dots, s ,
```

and include both endpoints ``\pm 1``.

The points are evaluated in the arithmetic `IT` and then converted to `T`. For a numeric `T`
the default is `IT=BigFloat`, so that the closed forms above are evaluated to full precision
whatever `T` is and the returned values are correctly rounded; see
[`ClenshawCurtisQuadrature`](@ref) for the trade-off involved in choosing a lower working
precision. For any other `T` the default is `IT=T`, so that with a symbolic `T` the closed
forms are evaluated exactly, ``\pi`` included; cf.
[`QuadratureRules._default_arithmetic`](@ref).

Points of the second kind require `s ≥ 2` and throw an `ErrorException` otherwise.

```jldoctest
julia> chebyshev_points(3, 1)
3-element Vector{Float64}:
 -0.8660254037844386
  0.0
  0.8660254037844386
```

See also [`chebyshev_nodes`](@ref) for the same points on ``[0,1]``.
"""
function chebyshev_points(::Type{T}, s::Integer, ::Val{1}; IT=_default_arithmetic(T)) where {T}
    T[ sin( IT(π) * (s-2i+1) / (2s) ) for i in s:-1:1 ]
end

function chebyshev_points(::Type{T}, s::Integer, ::Val{2}; IT=_default_arithmetic(T)) where {T}
    if s == 1
        throw(ErrorException("Chebyshev points of the second kind are not defined for one point."))
    end

    T[ cos( IT(π) * (i-1) / (s-1) ) for i in s:-1:1 ]
end

chebyshev_points(s, kind; kwargs...) = chebyshev_points(Float64, s, Val(kind); kwargs...)

@doc raw"""
    chebyshev_nodes(s, kind; IT=BigFloat)
    chebyshev_nodes(T, s, ::Val{kind}; IT=_default_arithmetic(T))

The `s` Chebyshev nodes of the first (`kind = 1`) or second (`kind = 2`) kind on the
interval ``[0,1]``, i.e., the Chebyshev points shifted and scaled from ``[-1,+1]``
to ``[0,1]``.

See [`chebyshev_points`](@ref) for the definition of the two kinds and for `IT`. Both the
points and the shift are computed in `IT` before the result is converted to `T`.

```jldoctest
julia> chebyshev_nodes(3, 1)
3-element Vector{Float64}:
 0.06698729810778067
 0.5
 0.9330127018922193
```
"""
function chebyshev_nodes(::Type{T}, s::Integer, kind::Val; IT=_default_arithmetic(T)) where {T}
    T.(shift_nodes(chebyshev_points(IT, s, kind; IT=IT)))
end

chebyshev_nodes(s, kind; kwargs...) = chebyshev_nodes(Float64, s, Val(kind); kwargs...)

@doc raw"""
    chebyshev_weights(s, kind; IT=BigFloat)
    chebyshev_weights(T, s, ::Val{kind}; IT=_default_arithmetic(T))

The `s` interpolatory weights for the interval ``[0,1]`` belonging to the Chebyshev nodes
of the first (`kind = 1`) or second (`kind = 2`) kind, summing to ``1``.

For `kind = 1` these are the weights of Fejér's first rule,

```math
w_i = \frac{1}{s} \left( 1 - 2 \sum_{j=1}^{\lfloor s/2 \rfloor}
      \frac{\cos (2 j \vartheta_i)}{4 j^2 - 1} \right) ,
      \qquad \vartheta_i = \frac{(2i-1) \pi}{2s} ,
```

already normalised to ``[0,1]``. For `kind = 2` the nodes coincide with the
Clenshaw-Curtis nodes, so the weights are those of
[`clenshaw_curtis_weights`](@ref), to which this function delegates.

See [`chebyshev_points`](@ref) for the definition of the two kinds and for the arguments.
`kind = 2` requires `s ≥ 2`.

```jldoctest
julia> chebyshev_weights(3, 1)
3-element Vector{Float64}:
 0.2222222222222222
 0.5555555555555556
 0.2222222222222222
```

See also [`chebyshev_point_weights`](@ref) for the same weights on ``[-1,+1]``.
"""
function chebyshev_weights(::Type{T}, s::Integer, ::Val{1}; IT=_default_arithmetic(T)) where {T}
    b = zeros(IT, s)
    for i in eachindex(b)
        tj = zero(IT)
        th = IT(π) * (2i-1) / (2s)
        for j in 1:div(s,2)
            tj += cos(2j*th) / IT(4j^2 - 1)
        end
        b[i] = (1 - 2tj) / IT(s)
    end

    T.(b)
end

chebyshev_weights(::Type{T}, s::Integer, ::Val{2}; kwargs...) where {T} = clenshaw_curtis_weights(T, s; kwargs...)

chebyshev_weights(s, kind; kwargs...) = chebyshev_weights(Float64, s, Val(kind); kwargs...)

@doc raw"""
    chebyshev_point_weights(s, kind; IT=BigFloat)
    chebyshev_point_weights(T, s, ::Val{kind}; IT=_default_arithmetic(T))

The `s` interpolatory weights for the interval ``[-1,+1]`` belonging to the Chebyshev
points of the first (`kind = 1`) or second (`kind = 2`) kind, i.e., the weights of
[`chebyshev_weights`](@ref) doubled so that they sum to ``2``.

See [`chebyshev_points`](@ref) for the definition of the two kinds and for the arguments.
`kind = 2` requires `s ≥ 2`.

```jldoctest
julia> chebyshev_point_weights(3, 1)
3-element Vector{Float64}:
 0.4444444444444444
 1.1111111111111112
 0.4444444444444444
```
"""
function chebyshev_point_weights(::Type{T}, s::Integer, kind::Val; IT=_default_arithmetic(T)) where {T}
    T.(unscale_weights(chebyshev_weights(IT, s, kind; IT=IT)))
end

chebyshev_point_weights(s, kind; kwargs...) = chebyshev_point_weights(Float64, s, Val(kind); kwargs...)


"""
    gauss_chebyshev_points(s; IT=BigFloat)
    gauss_chebyshev_points(T, s; IT=_default_arithmetic(T))

The `s` Gauss-Chebyshev points on the interval ``[-1,+1]``, i.e., the Chebyshev points
of the first kind, cf. [`chebyshev_points`](@ref).

These are the nodes of [`GaussChebyshevQuadrature`](@ref).
"""
gauss_chebyshev_points(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_points(T, s, Val(1); kwargs...)
gauss_chebyshev_points(s; kwargs...) = gauss_chebyshev_points(Float64, s; kwargs...)

"""
    gauss_chebyshev_nodes(s; IT=BigFloat)
    gauss_chebyshev_nodes(T, s; IT=_default_arithmetic(T))

The `s` Gauss-Chebyshev nodes on the interval ``[0,1]``, i.e., the Chebyshev nodes of
the first kind, cf. [`chebyshev_nodes`](@ref).
"""
gauss_chebyshev_nodes(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_nodes(T, s, Val(1); kwargs...)
gauss_chebyshev_nodes(s; kwargs...) = gauss_chebyshev_nodes(Float64, s; kwargs...)

"""
    gauss_chebyshev_weights(s; IT=BigFloat)
    gauss_chebyshev_weights(T, s; IT=_default_arithmetic(T))

The `s` Gauss-Chebyshev weights for the interval ``[0,1]``, i.e., the weights of Fejér's
first rule, cf. [`chebyshev_weights`](@ref). They sum to `1`.

These are the weights of [`GaussChebyshevQuadrature`](@ref).
"""
gauss_chebyshev_weights(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_weights(T, s, Val(1); kwargs...)
gauss_chebyshev_weights(s; kwargs...) = gauss_chebyshev_weights(Float64, s; kwargs...)

"""
    gauss_chebyshev_point_weights(s; IT=BigFloat)
    gauss_chebyshev_point_weights(T, s; IT=_default_arithmetic(T))

The `s` Gauss-Chebyshev weights for the interval ``[-1,+1]``, cf.
[`chebyshev_point_weights`](@ref). They sum to `2`.
"""
gauss_chebyshev_point_weights(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_point_weights(T, s, Val(1); kwargs...)
gauss_chebyshev_point_weights(s; kwargs...) = gauss_chebyshev_point_weights(Float64, s; kwargs...)

"""
    lobatto_chebyshev_points(s; IT=BigFloat)
    lobatto_chebyshev_points(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Chebyshev points on the interval ``[-1,+1]``, i.e., the Chebyshev points
of the second kind, cf. [`chebyshev_points`](@ref). They include both endpoints.

These points coincide with [`clenshaw_curtis_points`](@ref), and correspondingly
[`LobattoChebyshevQuadrature`](@ref) coincides with [`ClenshawCurtisQuadrature`](@ref).

Requires `s ≥ 2`.
"""
lobatto_chebyshev_points(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_points(T, s, Val(2); kwargs...)
lobatto_chebyshev_points(s; kwargs...) = lobatto_chebyshev_points(Float64, s; kwargs...)

"""
    lobatto_chebyshev_nodes(s; IT=BigFloat)
    lobatto_chebyshev_nodes(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Chebyshev nodes on the interval ``[0,1]``, i.e., the Chebyshev nodes of
the second kind, cf. [`chebyshev_nodes`](@ref). The first and last node are `0` and `1`.

Requires `s ≥ 2`.
"""
lobatto_chebyshev_nodes(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_nodes(T, s, Val(2); kwargs...)
lobatto_chebyshev_nodes(s; kwargs...) = lobatto_chebyshev_nodes(Float64, s; kwargs...)

"""
    lobatto_chebyshev_weights(s; IT=BigFloat)
    lobatto_chebyshev_weights(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Chebyshev weights for the interval ``[0,1]``, cf.
[`chebyshev_weights`](@ref). They sum to `1`.

As the nodes coincide with the Clenshaw-Curtis nodes, so do the weights: these are
identical to [`clenshaw_curtis_weights`](@ref), and they are the weights of
[`LobattoChebyshevQuadrature`](@ref).

Requires `s ≥ 2`.
"""
lobatto_chebyshev_weights(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_weights(T, s, Val(2); kwargs...)
lobatto_chebyshev_weights(s; kwargs...) = lobatto_chebyshev_weights(Float64, s; kwargs...)

"""
    lobatto_chebyshev_point_weights(s; IT=BigFloat)
    lobatto_chebyshev_point_weights(T, s; IT=_default_arithmetic(T))

The `s` Lobatto-Chebyshev weights for the interval ``[-1,+1]``, cf.
[`chebyshev_point_weights`](@ref). They sum to `2`.

Requires `s ≥ 2`.
"""
lobatto_chebyshev_point_weights(::Type{T}, s::Integer; kwargs...) where {T} = chebyshev_point_weights(T, s, Val(2); kwargs...)
lobatto_chebyshev_point_weights(s; kwargs...) = lobatto_chebyshev_point_weights(Float64, s; kwargs...)

@doc raw"""
    GaussChebyshevQuadrature(s; IT=BigFloat)
    GaussChebyshevQuadrature(T, s; IT=_default_arithmetic(T))

The Gauss-Chebyshev quadrature rule with `s` nodes on the interval ``[0,1]``, also known
as **Fejér's first rule**.

The nodes are the Chebyshev points of the first kind and the weights are those of the
interpolatory rule through them, obtained by integrating the Chebyshev interpolant of
`f` term by term:

```math
w_i = \frac{2}{s} \left( 1 - 2 \sum_{j=1}^{\lfloor s/2 \rfloor}
      \frac{\cos ( 2 j \theta_i )}{4 j^2 - 1} \right) ,
      \qquad \theta_i = \frac{(2i-1) \, \pi}{2s} ,
```

scaled by ``1/2`` for the interval ``[0,1]``. Being an interpolatory rule on `s` nodes,
it is exact for polynomials of degree ``\le s-1``. For odd `s` it gains one further
degree, because the monomial of degree `s` that it would otherwise fail on is odd about
the midpoint of the interval and its error therefore cancels. The **order is `s` for even
`s` and `s+1` for odd `s`**, and it is sharp. All weights are positive.

!!! note "This is not the weighted Gauss-Chebyshev rule"
    The classical Gauss-Chebyshev rule approximates the weighted integral
    ``\int_{-1}^{+1} f(x) \, (1-x^2)^{-1/2} \, dx`` with equal weights ``\pi/s`` and is
    exact to degree ``2s-1``. The rule implemented here uses the same nodes but
    approximates the *unweighted* integral ``\int_0^1 f(x) \, dx``, which is what
    [`QuadratureRule`](@ref) evaluates, and is therefore only exact to degree ``s-1``.

# Arguments
- `T`: element type of the resulting rule, `Float64` if omitted.
- `s`: number of nodes.
- `IT`: arithmetic in which nodes and weights are computed, `BigFloat` for a numeric `T` and
  `T` itself otherwise, cf. [`QuadratureRules._default_arithmetic`](@ref). As for
  [`ClenshawCurtisQuadrature`](@ref), the weight sum costs ``O(s^2)`` operations, so a lower
  working precision is faster but accumulates round-off in the intermediate terms; see the
  note there.

```jldoctest
julia> quad = GaussChebyshevQuadrature(3);

julia> order(quad)          # odd s, so exact up to degree 3
4

julia> quad(x -> x^2)
0.33333333333333337

julia> GaussChebyshevQuadrature(1) == MidpointQuadrature()
true
```

See also [`ClenshawCurtisQuadrature`](@ref) and [`ChebyshevQuadrature`](@ref).
"""
function GaussChebyshevQuadrature(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    c = chebyshev_nodes(IT, s, Val(1); IT=IT)
    b = chebyshev_weights(IT, s, Val(1); IT=IT)

    QuadratureRule(isodd(s) ? s+1 : s, c, b, T)
end

GaussChebyshevQuadrature(s; kwargs...) = GaussChebyshevQuadrature(Float64, s; kwargs...)


@doc raw"""
    LobattoChebyshevQuadrature(s; IT=BigFloat)
    LobattoChebyshevQuadrature(T, s; IT=_default_arithmetic(T))

The Lobatto-Chebyshev quadrature rule with `s` nodes on the interval ``[0,1]``.

Its nodes are the Chebyshev points of the second kind, which include both endpoints of
the interval. These are precisely the Clenshaw-Curtis nodes, and since an interpolatory
rule is uniquely determined by its nodes, this rule *is* the Clenshaw-Curtis rule:
`LobattoChebyshevQuadrature(s) == ClenshawCurtisQuadrature(s)`. The implementation
therefore delegates to [`ClenshawCurtisQuadrature`](@ref), where the weights and the
resulting **order**, `s` for even `s` and `s+1` for odd `s`, are documented.

Throws an `ErrorException` for `s == 1`.

```jldoctest
julia> LobattoChebyshevQuadrature(3) == ClenshawCurtisQuadrature(3)
true

julia> LobattoChebyshevQuadrature(3)
QuadratureRule{Float64, 3}(4, [0.0, 0.5, 1.0], [0.16666666666666666, 0.6666666666666666, 0.16666666666666666])
```
"""
function LobattoChebyshevQuadrature(::Type{T}, s::Integer; IT=_default_arithmetic(T)) where {T}
    if s == 1
        throw(ErrorException("Lobatto-Chebyshev quadrature is not defined for one stage."))
    end

    ClenshawCurtisQuadrature(T, s; IT=IT)
end

LobattoChebyshevQuadrature(s; kwargs...) = LobattoChebyshevQuadrature(Float64, s; kwargs...)


"""
    ChebyshevQuadrature(s, kind; kwargs...)
    ChebyshevQuadrature(T, s, ::Val{kind}; kwargs...)

Umbrella constructor for the Chebyshev quadrature rules with `s` nodes on ``[0,1]``,
dispatching on the kind of the Chebyshev points:

- `kind = 1` gives [`GaussChebyshevQuadrature`](@ref), i.e., Fejér's first rule,
- `kind = 2` gives [`LobattoChebyshevQuadrature`](@ref), i.e., the Clenshaw-Curtis rule.

Keyword arguments are forwarded to the selected rule; both accept `IT`, the arithmetic in
which nodes and weights are computed. An unsupported keyword raises a `MethodError` rather
than being silently ignored.

```jldoctest
julia> ChebyshevQuadrature(4, 1) == GaussChebyshevQuadrature(4)
true

julia> ChebyshevQuadrature(4, 2) == LobattoChebyshevQuadrature(4)
true

julia> ChebyshevQuadrature(8, 1; IT=Float64) ≈ ChebyshevQuadrature(8, 1)
true
```
"""
ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{1}; kwargs...) where {T} = GaussChebyshevQuadrature(T, s; kwargs...)
ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{2}; kwargs...) where {T} = LobattoChebyshevQuadrature(T, s; kwargs...)

ChebyshevQuadrature(s, kind; kwargs...) = ChebyshevQuadrature(Float64, s, Val(kind); kwargs...)
