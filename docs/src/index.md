```@meta
CurrentModule = QuadratureRules
```

# QuadratureRules


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
julia> TrapezoidalQuadrature()
QuadratureRule{Float64,2}(2, [0.0, 1.0], [0.5, 0.5])
```

The `QuadratureRule` type has the following fields:
- `order` the order of the method,
- `nodes` the nodes,
- `weights` the weights.

A functor is defined, which integrates some function `f(x)` using the quadrature rule:
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
- [`RiemannQuadratureLeft`](@ref)
- [`RiemannQuadratureRight`](@ref)
- [`MidpointQuadrature`](@ref)
- [`TrapezoidalQuadrature`](@ref)

as well as functions for generating quadrature rules on the fly:
- [`ClenshawCurtisQuadrature`](@ref)
- [`GaussChebyshevQuadrature`](@ref)
- [`GaussLegendreQuadrature`](@ref)
- [`LobattoChebyshevQuadrature`](@ref)
- [`LobattoLegendreQuadrature`](@ref)
