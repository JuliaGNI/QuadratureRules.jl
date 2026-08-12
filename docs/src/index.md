```@meta
CurrentModule = QuadratureRules
DocTestSetup = quote
    using QuadratureRules
end
```

# QuadratureRules

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
a [`QuadratureRule`](@ref) can be created by calling any one of the provided constructors, for example
```jldoctest intro
julia> quad = TrapezoidalQuadrature()
QuadratureRule{Float64, 2}(2, [0.0, 1.0], [0.5, 0.5])
```

The `QuadratureRule` type has the following fields:
- `order` the order of the method,
- `nodes` the nodes,
- `weights` the weights.

A functor is defined, which integrates functions `f(x)` over the interval ``[0,1]`` using the quadrature rule:
```jldoctest intro
julia> quad(x -> x^2)
0.5
```

There are several convenience functions for accessing the fields:
- [`nnodes`](@ref) the number of nodes,
- [`nodes`](@ref) the nodes,
- [`order`](@ref) the order,
- [`weights`](@ref) the weights,

as well as a function for looping over all nodes and weights:
- `Base.eachindex(quad::QuadratureRule) = eachindex(quad.nodes, quad.weights)`

All rules are defined on the reference interval ``[0,1]`` and their weights sum to one.
An integral over a general interval ``[a,b]`` is obtained by rescaling the argument:
```jldoctest intro
julia> a, b = 0.0, 2.0;

julia> quad = GaussLegendreQuadrature(8);

julia> (b - a) * quad(ξ -> exp(a + (b - a) * ξ))   # ∫₀² exp(x) dx = e² - 1
6.38905609893065

julia> exp(2) - 1
6.38905609893065
```


## Available Rules

There are several pre-tabulated quadrature rules:
- [`RiemannQuadratureLeft`](@ref)
- [`RiemannQuadratureRight`](@ref)
- [`MidpointQuadrature`](@ref)
- [`TrapezoidalQuadrature`](@ref)

as well as functions for generating quadrature rules with an arbitrary number of nodes on the fly:
- [`ChebyshevQuadrature`](@ref)
- [`ClenshawCurtisQuadrature`](@ref)
- [`GaussChebyshevQuadrature`](@ref)
- [`GaussLegendreQuadrature`](@ref)
- [`LobattoChebyshevQuadrature`](@ref)
- [`LobattoLegendreQuadrature`](@ref)
- [`TanhSinhQuadrature`](@ref)

Each generated rule accepts the number of nodes `s` and, optionally, the element type of
the resulting rule:
```jldoctest intro
julia> GaussLegendreQuadrature(2)
QuadratureRule{Float64, 2}(4, [0.2113248654051871, 0.7886751345948129], [0.5, 0.5])

julia> eltype(GaussLegendreQuadrature(BigFloat, 2))
BigFloat
```

Nodes and weights are computed in an internal working precision, controlled by the
keyword argument `IT` and defaulting to `BigFloat`, and converted to the requested
element type only at the end, so that the result is accurate to full precision. Where
this is not needed, the Legendre rules accept `fast=true`, which takes the nodes and
weights directly from
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
in double precision:
```jldoctest intro
julia> GaussLegendreQuadrature(5) ≈ GaussLegendreQuadrature(5; fast=true)
true
```

[`TanhSinhQuadrature`](@ref) is the exception to both statements. It is the trapezoidal rule
after a double-exponential change of variables, so its argument is a refinement level `n`
rather than a number of nodes — the latter follows from the level and the precision — and it
has no polynomial degree of exactness, its [`order`](@ref) being reported as `0`. What it
offers instead is an integrand singular at either endpoint:
```jldoctest intro
julia> TanhSinhQuadrature(3)(x -> log(x)) ≈ -1     # ∫₀¹ log x dx
true
```

See [Quadrature Rules](@ref) for a derivation of each rule.


## Accessing Nodes Directly

The nodes of each rule are also available without the corresponding weights. Functions
named `*_points` return them on the interval ``[-1,+1]``, functions named `*_nodes` on
the interval ``[0,1]``:
```jldoctest intro
julia> gauss_legendre_points(2)
2-element Vector{Float64}:
 -0.5773502691896257
  0.5773502691896257

julia> gauss_legendre_nodes(2)
2-element Vector{Float64}:
 0.2113248654051871
 0.7886751345948129
```

These exist for all rules, i.e., as [`gauss_legendre_points`](@ref),
[`lobatto_legendre_points`](@ref), [`chebyshev_points`](@ref),
[`gauss_chebyshev_points`](@ref), [`lobatto_chebyshev_points`](@ref),
[`clenshaw_curtis_points`](@ref) and [`tanh_sinh_points`](@ref), together with the
corresponding `*_nodes` functions.


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
