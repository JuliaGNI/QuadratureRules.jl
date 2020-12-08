# Quadrature Rules

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/QuadratureRules.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliagni.github.io/QuadratureRules.jl/dev)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/QuadratureRules.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/QuadratureRules.html)
[![Build Status](https://github.com/JuliaGNI/QuadratureRules.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/QuadratureRules.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/QuadratureRules.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/QuadratureRules.jl)
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
QuadratureRule{Float64,2}(2, [0.0, 1.0], [0.5, 0.5])
```

The `QuadratureRule` type has the following fields:
- `order` the order of the method,
- `nodes` the nodes,
- `weights` the weights.

A functor is defined, which integrates functions `f(x)` using the quadrature rule:
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


## Quadrature Rules

There are several pre-tabulated quadrature rules:
- `RiemannQuadratureLeft`
- `RiemannQuadratureRight`
- `MidpointQuadrature`
- `TrapezoidalQuadrature`

as well as functions for generating quadrature rules with an arbitrary number of nodes on the fly:
- `ClenshawCurtisQuadrature`
- `GaussChebyshevQuadrature`
- `GaussLegendreQuadrature`
- `LobattoChebyshevQuadrature`
- `LobattoLegendreQuadrature`


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
