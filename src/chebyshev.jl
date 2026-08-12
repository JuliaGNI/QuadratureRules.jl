
@doc raw"""
    chebyshev_points(s, kind)
    chebyshev_points(T, s, ::Val{kind})

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

and include both endpoints ``\pm 1``. All points are computed in `BigFloat` arithmetic
and converted to `T`.

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
function chebyshev_points(::Type{T}, s::Integer, ::Val{1}) where {T}
    T[ @big sin( π*(s-2i+1) / (2s) ) for i in s:-1:1 ]
end

function chebyshev_points(::Type{T}, s::Integer, ::Val{2}) where {T}
    if s == 1
        throw(ErrorException("Chebyshev points of the second kind are not defined for one point."))
    end

    T[ @big cos( π * (i-1) / (s-1) ) for i in s:-1:1 ]
end

chebyshev_points(s, kind) = chebyshev_points(Float64, s, Val(kind))

@doc raw"""
    chebyshev_nodes(s, kind)
    chebyshev_nodes(T, s, ::Val{kind})

The `s` Chebyshev nodes of the first (`kind = 1`) or second (`kind = 2`) kind on the
interval ``[0,1]``, i.e., the Chebyshev points shifted and scaled from ``[-1,+1]``
to ``[0,1]``.

See [`chebyshev_points`](@ref) for the definition of the two kinds.

```jldoctest
julia> chebyshev_nodes(3, 1)
3-element Vector{Float64}:
 0.0669872981077807
 0.5
 0.9330127018922193
```
"""
function chebyshev_nodes(::Type{T}, s::Integer, kind::Val) where {T}
    shift_nodes(chebyshev_points(T, s, kind))
end

chebyshev_nodes(s, kind) = chebyshev_nodes(Float64, s, Val(kind))


"""
    gauss_chebyshev_points(s)
    gauss_chebyshev_points(T, s)

The `s` Gauss-Chebyshev points on the interval ``[-1,+1]``, i.e., the Chebyshev points
of the first kind, cf. [`chebyshev_points`](@ref).

These are the nodes of [`GaussChebyshevQuadrature`](@ref).
"""
gauss_chebyshev_points(::Type{T}, s::Integer) where {T} = chebyshev_points(T, s, Val(1))
gauss_chebyshev_points(s) = gauss_chebyshev_points(Float64, s)

"""
    gauss_chebyshev_nodes(s)
    gauss_chebyshev_nodes(T, s)

The `s` Gauss-Chebyshev nodes on the interval ``[0,1]``, i.e., the Chebyshev nodes of
the first kind, cf. [`chebyshev_nodes`](@ref).
"""
gauss_chebyshev_nodes(::Type{T}, s::Integer) where {T} = chebyshev_nodes(T, s, Val(1))
gauss_chebyshev_nodes(s) = gauss_chebyshev_nodes(Float64, s)

"""
    lobatto_chebyshev_points(s)
    lobatto_chebyshev_points(T, s)

The `s` Lobatto-Chebyshev points on the interval ``[-1,+1]``, i.e., the Chebyshev points
of the second kind, cf. [`chebyshev_points`](@ref). They include both endpoints.

These points coincide with [`clenshaw_curtis_points`](@ref), and correspondingly
[`LobattoChebyshevQuadrature`](@ref) coincides with [`ClenshawCurtisQuadrature`](@ref).

Requires `s ≥ 2`.
"""
lobatto_chebyshev_points(::Type{T}, s::Integer) where {T} = chebyshev_points(T, s, Val(2))
lobatto_chebyshev_points(s) = lobatto_chebyshev_points(Float64, s)

"""
    lobatto_chebyshev_nodes(s)
    lobatto_chebyshev_nodes(T, s)

The `s` Lobatto-Chebyshev nodes on the interval ``[0,1]``, i.e., the Chebyshev nodes of
the second kind, cf. [`chebyshev_nodes`](@ref). The first and last node are `0` and `1`.

Requires `s ≥ 2`.
"""
lobatto_chebyshev_nodes(::Type{T}, s::Integer) where {T} = chebyshev_nodes(T, s, Val(2))
lobatto_chebyshev_nodes(s) = lobatto_chebyshev_nodes(Float64, s)

@doc raw"""
    GaussChebyshevQuadrature(s; IT=BigFloat)
    GaussChebyshevQuadrature(T, s; IT=BigFloat)

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
it is exact for polynomials of degree ``\le s-1``, so its **order is `s`**; for odd `s`
it gains one further degree by symmetry. All weights are positive.

!!! note "This is not the weighted Gauss-Chebyshev rule"
    The classical Gauss-Chebyshev rule approximates the weighted integral
    ``\int_{-1}^{+1} f(x) \, (1-x^2)^{-1/2} \, dx`` with equal weights ``\pi/s`` and is
    exact to degree ``2s-1``. The rule implemented here uses the same nodes but
    approximates the *unweighted* integral ``\int_0^1 f(x) \, dx``, which is what
    [`QuadratureRule`](@ref) evaluates, and is therefore only exact to degree ``s-1``.

# Arguments
- `T`: element type of the resulting rule, `Float64` if omitted.
- `s`: number of nodes.
- `IT`: arithmetic in which nodes and weights are computed, `BigFloat` by default. As for
  [`ClenshawCurtisQuadrature`](@ref), the weight sum costs ``O(s^2)`` operations, so
  `IT=Float64` is substantially faster when `T` is `Float64` anyway.

```jldoctest
julia> quad = GaussChebyshevQuadrature(3);

julia> order(quad)
3

julia> quad(x -> x^2)
0.33333333333333337
```

See also [`ClenshawCurtisQuadrature`](@ref) and [`ChebyshevQuadrature`](@ref).
"""
function GaussChebyshevQuadrature(::Type{T}, s::Integer; IT=BigFloat) where {T}
    c = chebyshev_nodes(IT, s, Val(1))
    b = zero(c)
    for i in eachindex(b)
        tj = zero(IT)
        th = IT(π) * (2i-1) / (2s)
        for j in 1:div(s,2)
            tj += cos(2j*th) / IT(4j^2 - 1)
        end
        b[i] = (1 - 2tj) / IT(s)
    end

    QuadratureRule(s, c, b, T)
end

GaussChebyshevQuadrature(s; kwargs...) = GaussChebyshevQuadrature(Float64, s; kwargs...)


@doc raw"""
    LobattoChebyshevQuadrature(s; IT=BigFloat)
    LobattoChebyshevQuadrature(T, s; IT=BigFloat)

The Lobatto-Chebyshev quadrature rule with `s` nodes on the interval ``[0,1]``.

Its nodes are the Chebyshev points of the second kind, which include both endpoints of
the interval. These are precisely the Clenshaw-Curtis nodes, and since an interpolatory
rule is uniquely determined by its nodes, this rule *is* the Clenshaw-Curtis rule:
`LobattoChebyshevQuadrature(s) == ClenshawCurtisQuadrature(s)`. The implementation
therefore delegates to [`ClenshawCurtisQuadrature`](@ref), where the weights and the
resulting **order `s`** are documented.

Throws an `ErrorException` for `s == 1`.

```jldoctest
julia> LobattoChebyshevQuadrature(3) == ClenshawCurtisQuadrature(3)
true

julia> LobattoChebyshevQuadrature(3)
QuadratureRule{Float64, 3}(3, [0.0, 0.5, 1.0], [0.16666666666666666, 0.6666666666666666, 0.16666666666666666])
```
"""
function LobattoChebyshevQuadrature(::Type{T}, s::Integer; IT=BigFloat) where {T}
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
