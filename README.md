# Quadrature Rules

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/QuadratureRules.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliagni.github.io/QuadratureRules.jl/latest)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/QuadratureRules.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/QuadratureRules.html)
[![Build Status](https://github.com/JuliaGNI/QuadratureRules.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/QuadratureRules.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/QuadratureRules.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/QuadratureRules.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.4310382.svg)](https://doi.org/10.5281/zenodo.4310382)


This package provides quadrature rules for numerical integration, e.g., in finite element methods or variational integrators. It provides a unified interface for quadrature rules from different sources and algorithms for the computation of quadrature rules with an arbitrary number of nodes and weights in arbitrary precision.

## Installation

*QuadratureRules.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```julia
]add QuadratureRules
```

## Basic Usage

After loading the Quadrature Rule module by
```julia
julia> using QuadratureRules
```
a `QuadratureRule` can be created by calling any one of the provided constructors, for example
```julia
julia> quad = TrapezoidalQuadrature()
QuadratureRule{Float64, 2}(2, [0.0, 1.0], [0.5, 0.5])
```

The `QuadratureRule` type has the following fields:
- `order` the order of the method,
- `nodes` the nodes,
- `weights` the weights.

A functor is defined, which integrates functions `f(x)` over the interval `[0,1]` using the quadrature rule:
```
julia> quad(x -> x^2)
0.5
```

There are several convenience functions for accessing the fields:
- `nnodes(::QuadratureRule{T,N}) where {T,N} = N`
- `nodes(quad::QuadratureRule) = quad.nodes`
- `order(quad::QuadratureRule) = quad.order`
- `weights(quad::QuadratureRule) = quad.weights`

as well as a function for looping over all nodes and weights:
- `Base.eachindex(quad::QuadratureRule) = eachindex(quad.nodes, quad.weights)`

All rules are defined on the reference interval `[0,1]` and their weights sum to one.
A rule of order `p` integrates polynomials of degree `≤ p-1` exactly.


## Quadrature Rules

There are several pre-tabulated quadrature rules:
- `RiemannQuadratureLeft`
- `RiemannQuadratureRight`
- `MidpointQuadrature`
- `TrapezoidalQuadrature`

as well as functions for generating quadrature rules with an arbitrary number of nodes on the fly:
- `ChebyshevQuadrature`
- `ClenshawCurtisQuadrature`
- `GaussChebyshevQuadrature`
- `GaussLegendreQuadrature`
- `LobattoChebyshevQuadrature`
- `LobattoLegendreQuadrature`
- `RadauLegendreQuadrature`
- `TanhSinhQuadrature`

Each of these takes the number of nodes `s` and, optionally, the element type of the
resulting rule, e.g. `GaussLegendreQuadrature(BigFloat, 5)`. Nodes and weights are
computed in an internal working precision, controlled by the keyword argument `IT` and
defaulting to `BigFloat`, and converted to the requested element type only at the end.
The Legendre rules additionally accept `fast=true`, which takes nodes and weights
directly from FastGaussQuadrature.jl in double precision.

`TanhSinhQuadrature` is the exception: it is the trapezoidal rule after a double-exponential
change of variables, so it takes a refinement level `n` rather than a number of nodes, and it
has no polynomial degree of exactness — its `order` is reported as `0`. Its nodes never reach
the endpoints, which makes it the rule to use for an integrand singular there, such as
`TanhSinhQuadrature(3)(log) ≈ -1`.


## Nodes and Weights

The nodes and weights of each rule are also available without constructing the rule. Each
family provides one node function and one weight function, both taking an `interval` keyword
argument that selects the interval the result lives on: `UnitInterval()`, the default, gives
`[0,1]`, where the weights sum to `1`, and `SymmetricInterval()` gives `[-1,+1]`, where they
sum to `2` — for every rule but tanh-sinh, which reaches those sums only up to its truncation
error:
```julia
julia> gauss_legendre_nodes(2)
2-element Vector{Float64}:
 0.2113248654051871
 0.7886751345948129

julia> gauss_legendre_weights(2)
2-element Vector{Float64}:
 0.5
 0.5

julia> gauss_legendre_nodes(2; interval = SymmetricInterval())
2-element Vector{Float64}:
 -0.5773502691896257
  0.5773502691896257
```

These exist for all rules, i.e., as `gauss_legendre_nodes`, `lobatto_legendre_nodes`,
`radau_legendre_nodes`, `chebyshev_nodes`, `gauss_chebyshev_nodes`,
`lobatto_chebyshev_nodes`, `clenshaw_curtis_nodes` and `tanh_sinh_nodes`, together with the
corresponding `*_weights` functions. Note that `SymmetricInterval` describes the interval, not
the node set — a Radau rule has deliberately asymmetric nodes on either interval.

The Radau rules take the prescribed endpoint as a further argument, `:left` or `:right`,
since the two variants are different rules:
```julia
julia> radau_legendre_nodes(2, :right)
2-element Vector{Float64}:
 0.3333333333333333
 1.0
```


## References

If you use QuadratureRules.jl in your work, please consider citing it by

```
@misc{Kraus:2020:QuadratureRules,
  title={QuadratureRules.jl: A Collection of Quadrature Rules in Julia},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/QuadratureRules.jl}},
  doi={10.5281/zenodo.4310382}
}
```


## Development

We are using git hooks, e.g., to enforce that all tests pass before pushing.
In order to activate these hooks, the following command must be executed once:
```
git config core.hooksPath .githooks
```
