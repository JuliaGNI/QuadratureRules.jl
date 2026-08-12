
@doc raw"""
    RiemannQuadratureLeft(T=Float64)

The left Riemann sum on ``[0,1]``, i.e., the one-node rule

```math
\int_0^1 f(x) \, dx \approx f(0) .
```

It is the interpolatory rule for the single node ``c_1 = 0`` and therefore exact for
constants only, giving order 1.

```jldoctest
julia> RiemannQuadratureLeft()
QuadratureRule{Float64, 1}(1, [0.0], [1.0])
```
"""
function RiemannQuadratureLeft(T=Float64)
    QuadratureRule(1, T[0], T[1])
end

@doc raw"""
    RiemannQuadratureRight(T=Float64)

The right Riemann sum on ``[0,1]``, i.e., the one-node rule

```math
\int_0^1 f(x) \, dx \approx f(1) .
```

It is the interpolatory rule for the single node ``c_1 = 1`` and therefore exact for
constants only, giving order 1.

```jldoctest
julia> RiemannQuadratureRight()
QuadratureRule{Float64, 1}(1, [1.0], [1.0])
```
"""
function RiemannQuadratureRight(T=Float64)
    QuadratureRule(1, T[1], T[1])
end

@doc raw"""
    MidpointQuadrature(T=Float64)

The midpoint rule on ``[0,1]``, i.e., the one-node rule

```math
\int_0^1 f(x) \, dx \approx f \big( \tfrac{1}{2} \big) .
```

Although it uses a single node, the symmetry of the node about the centre of the
interval makes it exact for linear functions as well, so its order is 2 rather than 1.
It is the one-node Gauss-Legendre rule, cf. [`GaussLegendreQuadrature`](@ref).

```jldoctest
julia> MidpointQuadrature()
QuadratureRule{Float64, 1}(2, [0.5], [1.0])
```
"""
function MidpointQuadrature(T=Float64)
    QuadratureRule(2, T[1//2], T[1])
end

@doc raw"""
    TrapezoidalQuadrature(T=Float64)

The trapezoidal rule on ``[0,1]``, i.e., the two-node rule

```math
\int_0^1 f(x) \, dx \approx \tfrac{1}{2} \, f(0) + \tfrac{1}{2} \, f(1) .
```

It integrates the linear interpolant through the endpoints and is thus exact for
linear functions, giving order 2. It is the two-node Lobatto-Legendre rule, cf.
[`LobattoLegendreQuadrature`](@ref).

```jldoctest
julia> TrapezoidalQuadrature()
QuadratureRule{Float64, 2}(2, [0.0, 1.0], [0.5, 0.5])
```
"""
function TrapezoidalQuadrature(T=Float64)
    QuadratureRule(2, T[0, 1], T[1//2, 1//2])
end
