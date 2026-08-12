# The Gauss-Legendre nodes on [-1,+1], i.e. the roots of P_s, in their own arithmetic. Taking
# the interval out of the way lets the accessor and the constructor share the root find.
function _gauss_legendre_nodes(s, IT)
    sort(_roots(_legendre_polynomial(s, IT), () -> FastGaussQuadrature.gausslegendre(s)[1]))
end

@doc raw"""
    gauss_legendre_nodes(s; kwargs...)
    gauss_legendre_nodes(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Gauss-Legendre nodes, i.e., the roots of the Legendre polynomial ``P_s`` mapped to
`interval`. They lie in the interior of the interval.

The roots are computed on ``[-1,+1]`` in the arithmetic `IT`, mapped to `interval` in that
same arithmetic and converted to `T` only at the very end. For a numeric `IT` they are
obtained by refining the double precision approximations of
`FastGaussQuadrature.gausslegendre` with Newton's method, so that `IT=BigFloat`, the default
whenever `T` is a numeric type, yields nodes that are accurate to full `BigFloat` precision
independently of `T`. For any other `IT`, in particular a symbolic one, they are instead
obtained exactly as the eigenvalues of the companion matrix of ``P_s``; this is the default
whenever `T` is not a numeric type, so that `gauss_legendre_nodes(T, s)` with a symbolic `T`
returns exact expressions.

# Arguments
- `T`: element type of the returned vector, `Float64` if omitted.
- `s`: number of nodes.
- `IT`: arithmetic in which the roots are computed, `BigFloat` for a numeric `T` and `T`
  itself otherwise, cf. [`QuadratureRules._default_arithmetic`](@ref).
- `interval`: [`UnitInterval`](@ref) for ``[0,1]``, the default and the interval of
  [`GaussLegendreQuadrature`](@ref), or [`SymmetricInterval`](@ref) for ``[-1,+1]``, related
  by ``c_i = (x_i + 1)/2``.

```jldoctest
julia> gauss_legendre_nodes(2)
2-element Vector{Float64}:
 0.2113248654051871
 0.7886751345948129

julia> gauss_legendre_nodes(2; interval = SymmetricInterval())
2-element Vector{Float64}:
 -0.5773502691896257
  0.5773502691896257
```
"""
function gauss_legendre_nodes(::Type{T}, s::Integer; IT=_default_arithmetic(T),
                              interval::QuadratureInterval=UnitInterval()) where {T}
    T.(_nodes_from_symmetric(_gauss_legendre_nodes(s, IT), interval))
end

gauss_legendre_nodes(s; kwargs...) = gauss_legendre_nodes(Float64, s; kwargs...)

# The Gauss-Legendre weights on [-1,+1] belonging to the precomputed nodes `x`, in their
# own arithmetic. Taking the nodes as an argument lets the constructor share the closed
# form with gauss_legendre_weights without repeating the root find.
function _gauss_legendre_weights(x::AbstractVector{IT}) where {IT}
    s = length(x)

    P = _legendre_polynomial(s, IT)
    D = Polynomials.derivative(P)

    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(IT[-x[i], 1]) )^2 )
        I(1) - I(-1)
    end

    [ inti(i) / D(x[i])^2  for i in 1:s ]
end

@doc raw"""
    gauss_legendre_weights(s; kwargs...)
    gauss_legendre_weights(T, s; IT=_default_arithmetic(T), interval=UnitInterval())

The `s` Gauss-Legendre weights belonging to the nodes returned by
[`gauss_legendre_nodes`](@ref) for the same `interval`. All of them are positive.

They are computed from the closed form

```math
w_i = \frac{1}{P_s'(x_i)^2} \int_{-1}^{+1} \left( \frac{P_s(x)}{x - x_i} \right)^2 dx ,
```

where the integrand is the square of the (unnormalised) Lagrange basis polynomial
associated with ``x_i``, evaluated by exact polynomial division and integration in the
arithmetic `IT`. Being formulated on ``[-1,+1]``, where the weights sum to ``2``, it is the
`SymmetricInterval` weights that are primary here; those on ``[0,1]`` are obtained as
``b_i = w_i / 2`` and sum to ``1``. See [`gauss_legendre_nodes`](@ref) for the arguments.

```jldoctest
julia> gauss_legendre_weights(2)
2-element Vector{Float64}:
 0.5
 0.5

julia> gauss_legendre_weights(2; interval = SymmetricInterval())
2-element Vector{Float64}:
 1.0
 1.0
```
"""
function gauss_legendre_weights(::Type{T}, s::Integer; IT=_default_arithmetic(T),
                                interval::QuadratureInterval=UnitInterval()) where {T}
    w = _gauss_legendre_weights(_gauss_legendre_nodes(s, IT))

    T.(_weights_from_symmetric(w, interval))
end

gauss_legendre_weights(s; kwargs...) = gauss_legendre_weights(Float64, s; kwargs...)


function _gauss_legendre_fast(s, T)
    c, b = FastGaussQuadrature.gausslegendre(s)
    shift!(b,c)
    QuadratureRule(2s, c, b, T)
end

@doc raw"""
    GaussLegendreQuadrature(s; IT=BigFloat, fast=false)
    GaussLegendreQuadrature(T, s; IT=_default_arithmetic(T), fast=false)

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
- `IT`: arithmetic in which nodes and weights are computed. For a numeric `T` this defaults
  to `BigFloat`, so that the result is accurate to full precision in `T`; for any other `T`,
  in particular a symbolic one, it defaults to `T` itself, so that the rule is computed
  exactly, cf. [`QuadratureRules._default_arithmetic`](@ref).
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
function GaussLegendreQuadrature(::Type{T}, s::Integer; IT=_default_arithmetic(T), fast=false) where {T}
    if fast
        return _gauss_legendre_fast(s, T)
    end

    x = _gauss_legendre_nodes(s, IT)
    w = _gauss_legendre_weights(x)

    return QuadratureRule(2s, shift_nodes(x), scale_weights(w), T)
end

GaussLegendreQuadrature(s; kwargs...) = GaussLegendreQuadrature(Float64, s; kwargs...)
