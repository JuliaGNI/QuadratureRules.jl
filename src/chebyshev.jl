# The Chebyshev nodes of the given kind on [-1,+1], in their own arithmetic. Both are
# generated in reverse index order so that the result is ascending.
function _chebyshev_nodes(s, ::Val{1}, IT)
    [sin(IT(π) * (s-2i+1) / (2s)) for i in s:-1:1]
end

function _chebyshev_nodes(s, ::Val{2}, IT)
    if s == 1
        throw(ErrorException("Chebyshev points of the second kind are not defined for one point."))
    end

    [cos(IT(π) * (i-1) / (s-1)) for i in s:-1:1]
end

@doc raw"""
    chebyshev_nodes(s, kind; kwargs...)
    chebyshev_nodes(T, s, ::Val{kind}; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Chebyshev nodes of the first (`kind = 1`) or second (`kind = 2`) kind, in ascending
order, mapped to `interval`.

On ``[-1,+1]`` the nodes of the first kind are the roots of the Chebyshev polynomial ``T_s``,

```math
x_i = \cos \left( \frac{(2i-1) \, \pi}{2s} \right) , \qquad i = 1, \dots, s ,
```

which lie strictly inside the interval. They are evaluated in the equivalent `sin` form
used in the implementation, which is more accurate near the ends of the interval. The
nodes of the second kind are the extrema of ``T_{s-1}``,

```math
x_i = \cos \left( \frac{(i-1) \, \pi}{s-1} \right) , \qquad i = 1, \dots, s ,
```

and include both endpoints.

The closed forms are evaluated on ``[-1,+1]`` in the arithmetic `IT`, mapped to `interval` in
that same arithmetic and converted to `T` only at the very end. For a numeric `T` the default
is `IT=BigFloat`, so that the closed forms are evaluated to full precision whatever `T` is and
the returned values are correctly rounded; see [`ClenshawCurtisQuadrature`](@ref) for the
trade-off involved in choosing a lower working precision. For any other `T` the default is
`IT=T`, so that with a symbolic `T` the closed forms are evaluated exactly, ``\pi`` included;
cf. [`QuadratureRules._default_arithmetic`](@ref).

Nodes of the second kind require `s ≥ 2` and throw an `ErrorException` otherwise.

# Arguments
- `T`: element type of the returned vector, `Float64` if omitted.
- `s`: number of nodes.
- `kind`: `1` or `2`, cf. above.
- `IT`: arithmetic in which the closed forms are evaluated, `BigFloat` for a numeric `T` and
  `T` itself otherwise, cf. [`QuadratureRules._default_arithmetic`](@ref).
- `interval`: [`UnitInterval`](@ref) for ``[0,1]``, the default and the interval of
  [`ChebyshevQuadrature`](@ref), or [`SymmetricInterval`](@ref) for ``[-1,+1]``.

```jldoctest
julia> chebyshev_nodes(3, 1)
3-element Vector{Float64}:
 0.06698729810778067
 0.5
 0.9330127018922193

julia> chebyshev_nodes(3, 1; interval = SymmetricInterval())
3-element Vector{Float64}:
 -0.8660254037844386
  0.0
  0.8660254037844386
```
"""
function chebyshev_nodes(::Type{T}, s::Integer, kind::Val; IT = _default_arithmetic(T),
        interval::QuadratureInterval = UnitInterval()) where {T}
    T.(_nodes_from_symmetric(_chebyshev_nodes(s, kind, IT), interval))
end

chebyshev_nodes(s, kind; kwargs...) = chebyshev_nodes(Float64, s, Val(kind); kwargs...)

# The Chebyshev weights of the given kind on [0,1], in their own arithmetic. Fejér's first
# rule for kind 1; for kind 2 the nodes are the Clenshaw-Curtis nodes, so the weights are too.
function _chebyshev_weights(s, ::Val{1}, IT)
    b = zeros(IT, s)
    for i in eachindex(b)
        tj = zero(IT)
        th = IT(π) * (2i-1) / (2s)
        for j in 1:div(s, 2)
            tj += cos(2j*th) / IT(4j^2 - 1)
        end
        b[i] = (1 - 2tj) / IT(s)
    end

    b
end

_chebyshev_weights(s, ::Val{2}, IT) = _clenshaw_curtis_weights(s, IT)

@doc raw"""
    chebyshev_weights(s, kind; kwargs...)
    chebyshev_weights(T, s, ::Val{kind}; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` interpolatory weights belonging to the Chebyshev nodes of the first (`kind = 1`) or
second (`kind = 2`) kind, i.e. to [`chebyshev_nodes`](@ref) for the same `kind` and
`interval`. All of them are positive.

For `kind = 1` these are the weights of Fejér's first rule,

```math
b_i = \frac{1}{s} \left( 1 - 2 \sum_{j=1}^{\lfloor s/2 \rfloor}
      \frac{\cos (2 j \vartheta_i)}{4 j^2 - 1} \right) ,
      \qquad \vartheta_i = \frac{(2i-1) \pi}{2s} ,
```

already normalised to ``[0,1]``, where they sum to ``1``. For `kind = 2` the nodes coincide
with the Clenshaw-Curtis nodes, so the weights are those of
[`clenshaw_curtis_weights`](@ref), to which this function delegates. Both formulas carry the
normalisation to ``[0,1]``, so unlike the Legendre families it is the `UnitInterval` weights
that are primary here; those on ``[-1,+1]`` are obtained as ``w_i = 2 b_i`` and sum to ``2``.

See [`chebyshev_nodes`](@ref) for the arguments. `kind = 2` requires `s ≥ 2`.

```jldoctest
julia> chebyshev_weights(3, 1)
3-element Vector{Float64}:
 0.2222222222222222
 0.5555555555555556
 0.2222222222222222

julia> chebyshev_weights(3, 1; interval = SymmetricInterval())
3-element Vector{Float64}:
 0.4444444444444444
 1.1111111111111112
 0.4444444444444444
```
"""
function chebyshev_weights(::Type{T}, s::Integer, kind::Val; IT = _default_arithmetic(T),
        interval::QuadratureInterval = UnitInterval()) where {T}
    T.(_weights_from_unit(_chebyshev_weights(s, kind, IT), interval))
end

chebyshev_weights(s, kind; kwargs...) = chebyshev_weights(Float64, s, Val(kind); kwargs...)

"""
    gauss_chebyshev_nodes(s; kwargs...)
    gauss_chebyshev_nodes(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Gauss-Chebyshev nodes, i.e., the Chebyshev nodes of the first kind, cf.
[`chebyshev_nodes`](@ref).

On the unit interval these are the nodes of [`GaussChebyshevQuadrature`](@ref).
"""
function gauss_chebyshev_nodes(::Type{T}, s::Integer; kwargs...) where {T}
    chebyshev_nodes(T, s, Val(1); kwargs...)
end
gauss_chebyshev_nodes(s; kwargs...) = gauss_chebyshev_nodes(Float64, s; kwargs...)

"""
    gauss_chebyshev_weights(s; kwargs...)
    gauss_chebyshev_weights(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Gauss-Chebyshev weights, i.e., the weights of Fejér's first rule, cf.
[`chebyshev_weights`](@ref).

On the unit interval these are the weights of [`GaussChebyshevQuadrature`](@ref).
"""
function gauss_chebyshev_weights(::Type{T}, s::Integer; kwargs...) where {T}
    chebyshev_weights(T, s, Val(1); kwargs...)
end
gauss_chebyshev_weights(s; kwargs...) = gauss_chebyshev_weights(Float64, s; kwargs...)

"""
    lobatto_chebyshev_nodes(s; kwargs...)
    lobatto_chebyshev_nodes(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Lobatto-Chebyshev nodes, i.e., the Chebyshev nodes of the second kind, cf.
[`chebyshev_nodes`](@ref). They include both endpoints of the interval.

These nodes coincide with [`clenshaw_curtis_nodes`](@ref), and correspondingly
[`LobattoChebyshevQuadrature`](@ref) coincides with [`ClenshawCurtisQuadrature`](@ref).

Requires `s ≥ 2`.
"""
function lobatto_chebyshev_nodes(::Type{T}, s::Integer; kwargs...) where {T}
    chebyshev_nodes(T, s, Val(2); kwargs...)
end
lobatto_chebyshev_nodes(s; kwargs...) = lobatto_chebyshev_nodes(Float64, s; kwargs...)

"""
    lobatto_chebyshev_weights(s; kwargs...)
    lobatto_chebyshev_weights(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Lobatto-Chebyshev weights, cf. [`chebyshev_weights`](@ref).

As the nodes coincide with the Clenshaw-Curtis nodes, so do the weights: these are
identical to [`clenshaw_curtis_weights`](@ref), and on the unit interval they are the
weights of [`LobattoChebyshevQuadrature`](@ref).

Requires `s ≥ 2`.
"""
function lobatto_chebyshev_weights(::Type{T}, s::Integer; kwargs...) where {T}
    chebyshev_weights(T, s, Val(2); kwargs...)
end
lobatto_chebyshev_weights(s; kwargs...) = lobatto_chebyshev_weights(Float64, s; kwargs...)

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
function GaussChebyshevQuadrature(::Type{T}, s::Integer; IT = _default_arithmetic(T)) where {T}
    c = _nodes_from_symmetric(_chebyshev_nodes(s, Val(1), IT), UnitInterval())
    b = _chebyshev_weights(s, Val(1), IT)

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
function LobattoChebyshevQuadrature(::Type{T}, s::Integer; IT = _default_arithmetic(T)) where {T}
    if s == 1
        throw(ErrorException("Lobatto-Chebyshev quadrature is not defined for one stage."))
    end

    ClenshawCurtisQuadrature(T, s; IT = IT)
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
function ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{1}; kwargs...) where {T}
    GaussChebyshevQuadrature(T, s; kwargs...)
end
function ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{2}; kwargs...) where {T}
    LobattoChebyshevQuadrature(T, s; kwargs...)
end

function ChebyshevQuadrature(s, kind; kwargs...)
    ChebyshevQuadrature(Float64, s, Val(kind); kwargs...)
end
