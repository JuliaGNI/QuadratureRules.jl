
"""
    QuadratureInterval

Supertype of the intervals on which nodes and weights can be requested, namely
[`UnitInterval`](@ref) and [`SymmetricInterval`](@ref).

Every node and weight function takes an `interval` keyword argument of this type, defaulting
to `UnitInterval()`. It selects the interval the result lives on and nothing else: the two
conventions describe the same quadrature rule, related by the affine map
``c_i = (x_i + 1)/2`` and ``b_i = w_i / 2``.
"""
abstract type QuadratureInterval end

"""
    UnitInterval()

The interval ``[0,1]``, on which the weights sum to ``1``.

This is the reference interval of [`QuadratureRule`](@ref) and the default for every node
and weight function. See [`QuadratureInterval`](@ref) and [`SymmetricInterval`](@ref).

```jldoctest
julia> gauss_legendre_nodes(2; interval = UnitInterval())
2-element Vector{Float64}:
 0.2113248654051871
 0.7886751345948129
```
"""
struct UnitInterval <: QuadratureInterval end

"""
    SymmetricInterval()

The interval ``[-1,+1]``, on which the weights sum to ``2``.

This is the interval on which the classical theory is formulated, being where the Legendre
and Chebyshev polynomials are defined, and the one the closed forms of most families are
computed on. It is symmetric about the origin; the nodes need not be, and for a Radau rule
deliberately are not. See [`QuadratureInterval`](@ref) and [`UnitInterval`](@ref).

```jldoctest
julia> gauss_legendre_nodes(2; interval = SymmetricInterval())
2-element Vector{Float64}:
 -0.5773502691896257
  0.5773502691896257
```
"""
struct SymmetricInterval <: QuadratureInterval end
