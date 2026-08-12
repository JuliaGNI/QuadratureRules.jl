
@doc raw"""
Nodes and weights of the tanh-sinh rule of level `n` on the interval ``[0,1]``, computed and
returned in the arithmetic `IT` and truncated to what the target type `T` can resolve. This
is the common back end of [`TanhSinhQuadrature`](@ref), [`tanh_sinh_nodes`](@ref),
[`tanh_sinh_points`](@ref), [`tanh_sinh_weights`](@ref) and
[`tanh_sinh_point_weights`](@ref), so that all five truncate identically.

Substituting ``x = \tanh ( \tfrac{\pi}{2} \sinh t )`` and applying the trapezoidal rule with
step ``h = 2^{-n}`` at ``t = k h`` gives, after the move to ``[0,1]``,

```math
\delta_k = \frac{1 - \tanh u_k}{2} = \frac{1}{1 + e^{2 u_k}} ,
\qquad
b_k = \frac{\pi h}{4} \, \frac{\cosh (k h)}{\cosh^2 u_k} ,
\qquad
u_k = \frac{\pi}{2} \sinh (k h) ,
```

where ``\delta_k`` is the distance of the outer pair of nodes from the endpoints,
``c_{\pm k} = \tfrac{1}{2} \pm (\tfrac{1}{2} - \delta_k)``. It is written through the
logistic form on the right to avoid the cancellation that ``1 - \tanh u_k`` would suffer for
large ``u_k`` — and it is those small numbers that the rule is all about.

The pairs are collected for ``k > 0`` and the rule is assembled symmetrically about the
centre node ``1/2``, so that the returned weights are symmetric bit for bit and the nodes
ascending by construction.

Level `n` uses ``h = 2^{-n}``, so the nodes of level `n` are contained in those of level
`n+1`. The single loop over ``k = 1, 2, \dots`` at the finest step is equivalent to refining
level by level: the truncation criterion is a threshold in ``t``, shared by all levels, so the
union of the level grids is exactly the multiples of ``2^{-n}`` below that threshold.
"""
function _tanh_sinh(::Type{T}, n::Integer; IT=BigFloat) where {T}
    if n < 1
        throw(ErrorException("Tanh-Sinh quadrature is not defined for less than one level."))
    end

    h = IT(1) / IT(2)^n

    d = IT[]
    b = IT[]

    k = 1

    while true
        t  = k * h
        u  = IT(π) / 2 * sinh(t)
        δₖ = inv(1 + exp(2u))
        bₖ = IT(π) / 4 * h * cosh(t) / cosh(u)^2

        # Truncate as soon as the pair can no longer be told apart from the endpoints, or
        # its weight from zero, in the target precision T. Both the node 1-δ on [0,1] and
        # the point 1-2δ on [-1,+1] are tested, so that neither representation degenerates.
        (iszero(T(bₖ)) || iszero(T(δₖ)) || isone(T(1 - δₖ)) || isone(T(1 - 2δₖ))) && break

        if !isempty(d) && (T(δₖ) == T(d[end]) || T(1 - δₖ) == T(1 - d[end]) ||
                           T(2δₖ - 1) == T(2 * d[end] - 1))
            # The pair is indistinguishable from its predecessor in T. Adding its weight to
            # that pair leaves the quadrature sum and the weight sum unchanged, and keeps
            # the rule strictly monotone in T. Doing it pairwise preserves the symmetry.
            b[end] += bₖ
        else
            push!(d, δₖ)
            push!(b, bₖ)
        end

        k += 1
    end

    c = [reverse(d); IT(1) / 2; 1 .- d]
    w = [reverse(b); IT(π) / 4 * h; b]

    return c, w
end


@doc raw"""
    tanh_sinh_points(n; IT=BigFloat)
    tanh_sinh_points(T, n; IT=BigFloat)

The tanh-sinh points of level `n` on the interval ``[-1,+1]``, in ascending order, i.e.
the images

```math
x_k = \tanh \left( \frac{\pi}{2} \sinh (k h) \right) ,
\qquad h = 2^{-n} ,
```

of the equidistant grid ``t = k h`` under the tanh-sinh transformation. They lie strictly
inside the interval and cluster double-exponentially fast towards its ends.

The number of points is *not* a free parameter: the grid is truncated where the points can
no longer be distinguished from ``\pm 1`` in the target precision, so it depends on `n`, on
`T` and, through the weights, on `IT`. See [`TanhSinhQuadrature`](@ref) for the arguments
and for the truncation criterion.

```jldoctest
julia> x = tanh_sinh_points(1)
13-element Vector{Float64}:
 -0.999999999999957
 -0.9999999888756649
 -0.9999774771924616
 -0.9975148564572244
 -0.9513679640727469
 -0.6742714922484359
  0.0
  0.6742714922484359
  0.9513679640727469
  0.9975148564572244
  0.9999774771924616
  0.9999999888756649
  0.999999999999957

julia> x == -reverse(x)
true
```

See also [`tanh_sinh_nodes`](@ref) for the same points on ``[0,1]``.
"""
function tanh_sinh_points(::Type{T}, n::Integer; IT=BigFloat) where {T}
    c, _ = _tanh_sinh(T, n; IT=IT)
    T.(unshift_nodes(c))
end

tanh_sinh_points(n; kwargs...) = tanh_sinh_points(Float64, n; kwargs...)

@doc raw"""
    tanh_sinh_nodes(n; IT=BigFloat)
    tanh_sinh_nodes(T, n; IT=BigFloat)

The tanh-sinh nodes of level `n` on the interval ``[0,1]``, i.e. the tanh-sinh points
shifted and scaled from ``[-1,+1]`` to ``[0,1]``. In this form the transformation is the
logistic sigmoid,

```math
c_k = \frac{1 + x_k}{2} = \frac{1}{1 + e^{-\pi \sinh (k h)}} ,
\qquad h = 2^{-n} ,
```

which is how the nodes are actually computed, since it yields the small distance of the
outer nodes from the endpoints without cancellation.

These are exactly the nodes of [`TanhSinhQuadrature`](@ref) at the same `n`, `T` and `IT`.
See [`tanh_sinh_points`](@ref) for the arguments.

```jldoctest
julia> tanh_sinh_nodes(1)
13-element Vector{Float64}:
 2.1470805279391204e-14
 5.562167559007666e-9
 1.1261403769203567e-5
 0.0012425717713878065
 0.024316017963626528
 0.1628642538757821
 0.5
 0.8371357461242179
 0.9756839820363735
 0.9987574282286122
 0.9999887385962308
 0.9999999944378325
 0.9999999999999786
```
"""
function tanh_sinh_nodes(::Type{T}, n::Integer; IT=BigFloat) where {T}
    c, _ = _tanh_sinh(T, n; IT=IT)
    T.(c)
end

tanh_sinh_nodes(n; kwargs...) = tanh_sinh_nodes(Float64, n; kwargs...)

@doc raw"""
    tanh_sinh_weights(n; IT=BigFloat)
    tanh_sinh_weights(T, n; IT=BigFloat)

The tanh-sinh weights of level `n` for the interval ``[0,1]``, belonging to the nodes
returned by [`tanh_sinh_nodes`](@ref), i.e. the trapezoidal weights ``h`` times the
Jacobian of the tanh-sinh transformation,

```math
b_k = \frac{\pi h}{4} \, \frac{\cosh (k h)}{\cosh^2 \big( \tfrac{\pi}{2} \sinh (k h) \big)} ,
\qquad h = 2^{-n} .
```

They are symmetric about the centre weight bit for bit, and all of them are positive.

Unlike every other family in this package these weights sum to ``1`` only up to the
truncation error, since the rule is not exact for the constant — see
[`TanhSinhQuadrature`](@ref), whose weights these are, for the arguments and for the
truncation criterion.

```jldoctest
julia> b = tanh_sinh_weights(1);

julia> b == reverse(b)
true

julia> isapprox(sum(b), 1; atol = 1E-5)
true
```

See also [`tanh_sinh_point_weights`](@ref) for the same weights on ``[-1,+1]``.
"""
function tanh_sinh_weights(::Type{T}, n::Integer; IT=BigFloat) where {T}
    _, b = _tanh_sinh(T, n; IT=IT)
    T.(b)
end

tanh_sinh_weights(n; kwargs...) = tanh_sinh_weights(Float64, n; kwargs...)

@doc raw"""
    tanh_sinh_point_weights(n; IT=BigFloat)
    tanh_sinh_point_weights(T, n; IT=BigFloat)

The tanh-sinh weights of level `n` for the interval ``[-1,+1]``, i.e. the weights of
[`tanh_sinh_weights`](@ref) doubled, belonging to the points returned by
[`tanh_sinh_points`](@ref). They sum to ``2`` only up to the truncation error.

See [`TanhSinhQuadrature`](@ref) for the arguments.

```jldoctest
julia> isapprox(sum(tanh_sinh_point_weights(1)), 2; atol = 1E-5)
true
```
"""
function tanh_sinh_point_weights(::Type{T}, n::Integer; IT=BigFloat) where {T}
    _, b = _tanh_sinh(T, n; IT=IT)
    T.(unscale_weights(b))
end

tanh_sinh_point_weights(n; kwargs...) = tanh_sinh_point_weights(Float64, n; kwargs...)


@doc raw"""
    TanhSinhQuadrature(n; IT=BigFloat)
    TanhSinhQuadrature(T, n; IT=BigFloat)

The tanh-sinh quadrature rule of level `n` on the interval ``[0,1]``.

Also known as the *double-exponential* formula of [takahasi1974](@citet). It is not an
interpolatory rule on a prescribed set of nodes, but the trapezoidal rule applied after the
change of variables

```math
x = \tanh \left( \frac{\pi}{2} \sinh t \right) ,
\qquad t \in \mathbb{R} ,
```

which maps ``\mathbb{R}`` onto ``(-1,+1)``. Since ``dx/dt`` decays like
``\exp ( - \tfrac{\pi}{2} e^{|t|} )``, the transformed integrand vanishes
double-exponentially at both ends whatever the integrand does at the endpoints, and the
infinite trapezoidal sum can be truncated after a handful of terms. With ``h = 2^{-n}`` the
nodes and weights on ``[0,1]`` are

```math
c_k = \frac{1}{1 + e^{-\pi \sinh (k h)}} ,
\qquad
b_k = \frac{\pi h}{4} \, \frac{\cosh (k h)}{\cosh^2 \big( \tfrac{\pi}{2} \sinh (k h) \big)} ,
\qquad k \in \mathbb{Z} ,
```

symmetric about the centre node ``c_0 = 1/2``. See the
[Tanh-Sinh](@ref "Tanh-Sinh quadrature") section of the manual for the derivation and for
the literature.

The rule is exact for **no** polynomial, not even for a constant: the weights sum to one
only up to the truncation error. Its [`order`](@ref) is therefore reported as `0`, and it is
the one rule in this package whose accuracy is not described by a polynomial degree of
exactness. What it does instead is converge like ``\mathcal{O} ( e^{-cN/\log N} )`` in the
number of nodes ``N``, which in practice means that **each level roughly doubles the number
of correct digits** until the precision of `T` is reached. For `Float64` that point is
reached at level 3; a `BigFloat` rule at the default precision needs level 5.

# Arguments
- `T`: element type of the resulting rule, `Float64` if omitted.
- `n`: the level, i.e. the number of halvings of the step size, so that ``h = 2^{-n}``.
  This is **not** a node count — the number of nodes follows from the truncation below and
  grows like ``2^n``, e.g. 13, 25, 51, 101, 203 for levels 1 to 5 at `T=Float64`.
- `IT`: arithmetic in which nodes and weights are computed, `BigFloat` by default.

The grid is truncated at the first `k` whose weight rounds to zero in `T` or whose node
rounds to an endpoint — of either ``[0,1]`` or ``[-1,+1]``, so that neither representation
degenerates. A pair that is indistinguishable from its predecessor in `T` has its weight
folded into that predecessor instead of being added, which leaves the quadrature sum
unchanged and keeps the nodes strictly increasing. Consequently no node ever coincides with
an endpoint, and integrands that are singular there may be passed in directly.

Throws an `ErrorException` for `n < 1`.

!!! note "Endpoint singularities and the precision of `T`"
    Tanh-sinh is the method of choice for an integrand with a singularity at an endpoint,
    but how far it can get is limited by how closely a node of type `T` can approach that
    endpoint, namely to within about `eps(T)`. For an integrand behaving like
    ``x^{-1/2}`` the neglected tail is therefore of size ``\sqrt{\texttt{eps(T)}}``, and no
    level beyond the third improves on that: about 8 correct digits in `Float64`, 39 at the
    default `BigFloat` precision, 78 at twice that. A logarithmic singularity is far
    milder and is integrated to full precision. This is the clearest illustration of why
    this package computes in arbitrary precision.

```jldoctest
julia> quad = TanhSinhQuadrature(3);

julia> nnodes(quad)
51

julia> order(quad)
0

julia> quad(x -> log(x)) ≈ -1                     # ∫₀¹ log x dx, singular at x = 0
true

julia> quad(x -> exp(x)) ≈ exp(1) - 1
true

julia> isfinite(quad(x -> 1 / sqrt(x * (1 - x))))  # no node sits on an endpoint
true
```

See also [`tanh_sinh_nodes`](@ref) and [`tanh_sinh_points`](@ref) for the nodes alone,
[`tanh_sinh_weights`](@ref) and [`tanh_sinh_point_weights`](@ref) for the weights, and
[`GaussLegendreQuadrature`](@ref), which is the better choice for an integrand that is
smooth up to and including the endpoints.
"""
function TanhSinhQuadrature(::Type{T}, n::Integer; IT=BigFloat) where {T}
    c, b = _tanh_sinh(T, n; IT=IT)

    QuadratureRule(0, c, b, T)
end

TanhSinhQuadrature(n; kwargs...) = TanhSinhQuadrature(Float64, n; kwargs...)
