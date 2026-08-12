```@meta
CurrentModule = QuadratureRules
```

# Library

The complete API of *QuadratureRules.jl*. See [Quadrature Rules](@ref) for how the
individual rules are derived and [Numerical Quadrature](@ref) for the underlying theory.


## The QuadratureRule type

Every rule provided by the package is an instance of the single concrete type
`QuadratureRule`. The rules themselves are ordinary functions returning such an
instance, not distinct types.

The following covers the type itself, its outer constructor, and the functor that
evaluates the quadrature sum for a given integrand.

```@docs
QuadratureRule
```


## Generic API

These functions grant access to the fields of a rule and implement the usual comparison
and iteration protocols.

```@docs
nnodes
nodes
order
weights
```

### Base methods

```@docs
Base.eachindex(::QuadratureRule)
Base.eltype(::QuadratureRule)
Base.:(==)(::QuadratureRule, ::QuadratureRule)
Base.isequal(::QuadratureRule, ::QuadratureRule)
Base.isapprox(::QuadratureRule, ::QuadratureRule)
```


## Tabulated quadrature rules

Classical low-order rules with a fixed number of nodes. Each takes an optional element
type as its only argument.

```@docs
RiemannQuadratureLeft
RiemannQuadratureRight
MidpointQuadrature
TrapezoidalQuadrature
```


## Generated quadrature rules

Rules computed on the fly for an arbitrary number of nodes `s` and in arbitrary
precision. Each has a method taking an explicit element type `T` and one defaulting to
`Float64`.

```@docs
GaussLegendreQuadrature
LobattoLegendreQuadrature
RadauLegendreQuadrature
ClenshawCurtisQuadrature
GaussChebyshevQuadrature
LobattoChebyshevQuadrature
ChebyshevQuadrature
```

[`TanhSinhQuadrature`](@ref) belongs here too, but differs from the rules above in two
respects: its argument is a refinement level `n` rather than a number of nodes, and it has no
polynomial degree of exactness, so its [`order`](@ref) is `0`.

```@docs
TanhSinhQuadrature
```


## Points, Nodes and Weights

The nodes and weights of every rule are also available on their own, without having to
construct the rule. Two conventions are offered throughout:

- the `*_points` and `*_point_weights` functions work on the interval ``[-1,+1]``, the
  interval on which the classical theory is formulated, and their weights sum to ``2``,
- the `*_nodes` and `*_weights` functions work on the interval ``[0,1]``, the reference
  interval used by [`QuadratureRule`](@ref), and their weights sum to ``1``,

related by ``c_i = (x_i + 1)/2`` and ``b_i = w_i / 2``. Each function has a method taking
an explicit element type `T` and one defaulting to `Float64`.

Within each family one convention is primary and the other is derived from it, so that the
two agree exactly at equal working precision and up to rounding across precisions. The
nodes are always derived from the points. For the weights it depends on the family: the
Legendre closed forms are formulated on ``[-1,+1]``, so there the `*_point_weights` are
primary, whereas the Chebyshev, Clenshaw-Curtis and tanh-sinh formulas already carry the
normalisation to ``[0,1]``, so there the `*_weights` are.

The weight sums quoted above hold for every family except tanh-sinh, which is exact for no
polynomial at all and therefore attains them only up to its truncation error.

### Legendre

```@docs
gauss_legendre_points
gauss_legendre_nodes
gauss_legendre_point_weights
gauss_legendre_weights
lobatto_legendre_points
lobatto_legendre_nodes
lobatto_legendre_point_weights
lobatto_legendre_weights
radau_legendre_points
radau_legendre_nodes
radau_legendre_point_weights
radau_legendre_weights
```

### Chebyshev

```@docs
chebyshev_points
chebyshev_nodes
chebyshev_point_weights
chebyshev_weights
gauss_chebyshev_points
gauss_chebyshev_nodes
gauss_chebyshev_point_weights
gauss_chebyshev_weights
lobatto_chebyshev_points
lobatto_chebyshev_nodes
lobatto_chebyshev_point_weights
lobatto_chebyshev_weights
clenshaw_curtis_points
clenshaw_curtis_nodes
clenshaw_curtis_point_weights
clenshaw_curtis_weights
```

### Tanh-Sinh

```@docs
tanh_sinh_points
tanh_sinh_nodes
tanh_sinh_point_weights
tanh_sinh_weights
```


## Internals

The following are not exported and are documented for reference only. They are not part
of the public API and may change without notice.

```@docs
QuadratureRules._default_arithmetic
QuadratureRules._legendre
QuadratureRules._legendre_polynomial
QuadratureRules._newton_roots
QuadratureRules._roots
QuadratureRules._tanh_sinh
QuadratureRules.shift_nodes
QuadratureRules.unshift_nodes
QuadratureRules.scale_weights
QuadratureRules.unscale_weights
QuadratureRules.shift!
QuadratureRules.unshift!
```
