
@doc raw"""
    QuadratureRule{T,N}

A quadrature rule on the reference interval ``[0,1]``, that is, an approximation of the form

```math
\int_0^1 f(x) \, dx \approx \sum_{i=1}^{N} b_i \, f(c_i) ,
```

with nodes ``c_i`` and weights ``b_i``. All rules in this package are normalised to
``[0,1]``, so that the weights sum to one. An integral over a general interval
``[a,b]`` is obtained by the affine change of variables

```math
\int_a^b f(x) \, dx = (b-a) \int_0^1 f \big( a + (b-a) \, \xi \big) \, d\xi .
```

# Type parameters
- `T`: element type of the nodes and weights, e.g. `Float64` or `BigFloat`.
- `N`: number of nodes, encoded in the type so that it is available at compile time.

# Fields
- `order::Int`: the order ``p`` of the rule; a rule of order ``p`` integrates polynomials
  of degree ``\le p-1`` exactly (see [`order`](@ref)).
- `nodes::Vector{T}`: the nodes ``c_i \in [0,1]``.
- `weights::Vector{T}`: the weights ``b_i``.

# Examples

A rule is usually obtained from one of the constructors listed in the
[Quadrature Rules](@ref) section rather than built by hand:

```jldoctest
julia> quad = TrapezoidalQuadrature()
QuadratureRule{Float64, 2}(2, [0.0, 1.0], [0.5, 0.5])

julia> quad(x -> x^2)
0.5
```

See also [`nodes`](@ref), [`weights`](@ref), [`order`](@ref), [`nnodes`](@ref).
"""
struct QuadratureRule{T,N}
    order::Int
    nodes::Vector{T}
    weights::Vector{T}
end

"""
    QuadratureRule(order, nodes, weights, TT=eltype(nodes))

Construct a [`QuadratureRule`](@ref) from `order`, `nodes` and `weights`.

The optional fourth positional argument `TT` is the element type of the resulting rule.
It allows the nodes and weights to be computed in a high-precision arithmetic and stored
in a lower-precision one, which is how the generated rules obtain accurately rounded
`Float64` values:

```jldoctest
julia> QuadratureRule(2, BigFloat[0, 1], BigFloat[1//2, 1//2], Float64)
QuadratureRule{Float64, 2}(2, [0.0, 1.0], [0.5, 0.5])
```

`nodes` and `weights` must have the same length and the same element type.
"""
function QuadratureRule(order, nodes::AbstractVector{T}, weights::AbstractVector{T}, TT=T) where {T}
    @assert length(nodes) == length(weights)
    QuadratureRule{TT, length(nodes)}(order, nodes, weights)
end

@doc raw"""
    (quad::QuadratureRule)(f)

Integrate the function `f` over the interval ``[0,1]``, i.e., evaluate
``\sum_i b_i \, f(c_i) \approx \int_0^1 f(x) \, dx``.

```jldoctest
julia> quad = GaussLegendreQuadrature(2);

julia> quad(x -> x^3)
0.25
```

To integrate over a general interval `[a,b]`, rescale the argument and multiply by the
length of the interval:

```jldoctest
julia> a, b = 1.0, 3.0;

julia> quad = GaussLegendreQuadrature(4);

julia> (b - a) * quad(ξ -> (a + (b - a) * ξ)^2)   # ∫₁³ x² dx = 26/3
8.666666666666668
```
"""
function (quad::QuadratureRule)(f::Function)
    sum(w * f(x) for (w, x) in zip(quad.weights, quad.nodes))
end


"""
    nnodes(quad::QuadratureRule)

Return the number of nodes of `quad`.

This is the type parameter `N` of [`QuadratureRule`](@ref) and therefore known at
compile time.

```jldoctest
julia> nnodes(GaussLegendreQuadrature(3))
3
```
"""
nnodes(::QuadratureRule{T,N}) where {T,N} = N

@doc raw"""
    nodes(quad::QuadratureRule)

Return the vector of nodes ``c_i \in [0,1]`` of `quad`.

```jldoctest
julia> nodes(TrapezoidalQuadrature())
2-element Vector{Float64}:
 0.0
 1.0
```
"""
nodes(quad::QuadratureRule) = quad.nodes

@doc raw"""
    order(quad::QuadratureRule)

Return the order ``p`` of `quad`.

A quadrature rule has order ``p`` if it integrates all polynomials of degree
``\le p-1`` exactly, that is, if the moment conditions

```math
\sum_{i} b_i \, c_i^{k} = \frac{1}{k+1} , \qquad k = 0, 1, \dots, p-1 ,
```

are satisfied. The order is a *guaranteed* value: the Clenshaw-Curtis and Chebyshev
rules gain one additional degree when the number of nodes is odd, which is not
reflected here.

```jldoctest
julia> order(GaussLegendreQuadrature(3))     # exact up to degree 5
6

julia> order(LobattoLegendreQuadrature(3))   # exact up to degree 3
4
```
"""
order(quad::QuadratureRule) = quad.order

@doc raw"""
    weights(quad::QuadratureRule)

Return the vector of weights ``b_i`` of `quad`.

As all rules are normalised to the interval ``[0,1]``, the weights sum to one.

```jldoctest
julia> weights(TrapezoidalQuadrature())
2-element Vector{Float64}:
 0.5
 0.5
```
"""
weights(quad::QuadratureRule) = quad.weights

"""
    eachindex(quad::QuadratureRule)

Return an iterator over the indices of the nodes and weights of `quad`, so that
all pairs ``(b_i, c_i)`` can be visited in a loop.

```jldoctest
julia> quad = TrapezoidalQuadrature();

julia> for i in eachindex(quad)
           println(nodes(quad)[i], " ", weights(quad)[i])
       end
0.0 0.5
1.0 0.5
```
"""
Base.eachindex(quad::QuadratureRule) = eachindex(quad.nodes, quad.weights)

Base.hash(quad::QuadratureRule, h::UInt) = hash(quad.order, hash(quad.nodes, hash(quad.weights, hash(:QuadratureRule, h))))

"""
    ==(quad1::QuadratureRule, quad2::QuadratureRule)

Compare two quadrature rules by value, i.e., check that their orders, nodes and weights
agree. Rules of different element type may compare equal; use [`isequal`](@ref) to
distinguish them.
"""
Base.:(==)(quad1::QuadratureRule, quad2::QuadratureRule) = (
                   quad1.order == quad2.order
                && quad1.nodes == quad2.nodes
                && quad1.weights == quad2.weights)

"""
    isequal(quad1::QuadratureRule, quad2::QuadratureRule)

Check that two quadrature rules are equal by value *and* of the same type, that is,
that they agree in their element type `T` and their number of nodes `N` in addition
to their orders, nodes and weights.
"""
Base.isequal(quad1::QuadratureRule{T1,N1}, quad2::QuadratureRule{T2,N2}) where {T1,N1,T2,N2} = (
                quad1 == quad2
                && T1 == T2
                && N1 == N2)

"""
    isapprox(quad1::QuadratureRule, quad2::QuadratureRule; kwargs...)

Check that two quadrature rules have the same order and approximately equal nodes and
weights. All keyword arguments are forwarded to `isapprox` for the nodes and weights.

This is the appropriate comparison for rules computed by different algorithms or in
different working precisions, e.g. with and without `fast=true`.
"""
Base.isapprox(quad1::QuadratureRule, quad2::QuadratureRule; kwargs...) = (
                            quad1.order == quad2.order
                && isapprox(quad1.nodes,   quad2.nodes;   kwargs...)
                && isapprox(quad1.weights, quad2.weights; kwargs...))

"""
    eltype(quad::QuadratureRule)

Return the element type `T` of the nodes and weights of `quad`.

```jldoctest
julia> eltype(GaussLegendreQuadrature(BigFloat, 3))
BigFloat
```
"""
Base.eltype(::QuadratureRule{T}) where {T} = T
