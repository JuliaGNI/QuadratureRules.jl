```@meta
CurrentModule = QuadratureRules
```

# Quadrature Rules

This page describes how each of the rules provided by the package is derived and
computed. All rules live on the reference interval $[0,1]$ and are returned as a
[`QuadratureRule`](@ref); see [Numerical Quadrature](@ref) for the underlying theory and
for the meaning of the reported order.

Two groups can be distinguished. The *tabulated* rules have a fixed, small number of
nodes and are written out explicitly. The *generated* rules are computed on the fly for
an arbitrary number of nodes `s` and in arbitrary precision.


## Tabulated rules

These are the classical low-order rules. Each is the interpolatory rule for its node set,
so its weights are the integrals of the corresponding Lagrange basis functions over
$[0,1]$.

### Riemann sums

Sampling at a single endpoint and integrating the resulting constant gives the left and
right Riemann sums,

```math
\int_0^1 f(x) \, dx \approx f(0) ,
\qquad
\int_0^1 f(x) \, dx \approx f(1) .
```

Both reproduce constants and nothing more, so their order is 1.

```@repl rules
using QuadratureRules
RiemannQuadratureLeft()
RiemannQuadratureRight()
```

### Midpoint rule

Sampling at the centre of the interval instead,

```math
\int_0^1 f(x) \, dx \approx f \big( \tfrac{1}{2} \big) ,
```

also integrates the constant interpolant, but the symmetry of the node about the centre
buys an extra degree: the error term for a linear function is odd about $x = 1/2$ and
integrates to zero. The midpoint rule is therefore exact for linear functions and has
order 2, despite using only one node. It is in fact the one-node Gauss-Legendre rule.

```@repl rules
MidpointQuadrature() == GaussLegendreQuadrature(1)
```

### Trapezoidal rule

Interpolating linearly between the two endpoints and integrating gives

```math
\int_0^1 f(x) \, dx \approx \tfrac{1}{2} \, f(0) + \tfrac{1}{2} \, f(1) ,
```

which is exact for linear functions and thus of order 2. It is the two-node
Lobatto-Legendre rule.

```@repl rules
TrapezoidalQuadrature() == LobattoLegendreQuadrature(2)
```


## Gauss-Legendre quadrature

[`GaussLegendreQuadrature`](@ref) is the rule of maximal degree of exactness: with `s`
nodes it integrates polynomials up to degree $2s-1$ exactly, so its order is $2s$.

### Nodes

The nodes are the roots of the Legendre polynomial $P_s$, mapped from $[-1,+1]$ to
$[0,1]$. The polynomial itself is built from the three-term recurrence

```math
j \, P_j (x) = (2j-1) \, x \, P_{j-1} (x) - (j-1) \, P_{j-2} (x) ,
\qquad P_0 = 1 , \quad P_1 = x ,
```

which is implemented by [`QuadratureRules._legendre`](@ref). The same routine is used
both to evaluate $P_j$ at a number and to construct it as a `Polynomial`, by starting the
recurrence from the polynomial $x$ instead of a scalar.

Its roots have no closed form. They are obtained by taking the double precision values
from `FastGaussQuadrature.gausslegendre` — which for small `s` computes them by the
Golub-Welsch eigenvalue algorithm [golub1969](@cite) — as initial guesses, and refining
them with Newton's method in the working precision `IT`
([`QuadratureRules._newton_roots`](@ref)), which iterates until the correction stops
decreasing. This yields nodes accurate to full `BigFloat` precision.

### Weights

Rather than integrating the Lagrange basis functions directly, the implementation uses
the closed form

```math
w_i = \frac{1}{\big[ P_s' (x_i) \big]^2}
      \int_{-1}^{+1} \left( \frac{P_s (x)}{x - x_i} \right)^{2} dx .
```

The quotient $P_s(x) / (x - x_i)$ is exactly the unnormalised Lagrange basis polynomial
belonging to $x_i$, since $P_s$ vanishes at all nodes; dividing by $P_s'(x_i)$ normalises
it to one at $x_i$. The division is carried out exactly with `Polynomials.÷`, the square
is integrated symbolically, and the result is evaluated at the endpoints — all in the
arithmetic `IT`, so no accuracy is lost. Finally the nodes are shifted and the weights
halved to move the rule to $[0,1]$.

```@repl rules
quad = GaussLegendreQuadrature(3)
quad(x -> x^5)     # exact up to degree 5
```

Passing `fast=true` bypasses this construction and takes both nodes and weights straight
from FastGaussQuadrature.jl in double precision. The result is accurate to about
`Float64` precision and agrees with the default path to that tolerance.

```@repl rules
GaussLegendreQuadrature(5) ≈ GaussLegendreQuadrature(5; fast=true)
```


## Lobatto-Legendre quadrature

[`LobattoLegendreQuadrature`](@ref) constrains both endpoints of the interval to be
nodes. Only $s-2$ nodes remain free, so the degree of exactness drops to $2s-3$ and the
order to $2s-2$. At least two nodes are required; `s == 1` throws an `ErrorException`.

### Nodes

The interior nodes are the roots of $P_{s-1}'$, together with the endpoints $\pm 1$.
Instead of differentiating the Legendre polynomial, the implementation exploits
Rodrigues' formula: $P_{s-1}$ is proportional to the $(s-1)$-st derivative of
$(x^2-1)^{s-1}$, so $P_{s-1}'$ has the same roots as the $(s-2)$-nd derivative of
$(1-x^2)^{s-1}$. That derivative is formed directly with `Polynomials.derivative`, which
avoids building and differentiating the Legendre polynomial itself.

The roots are again Newton-refined from the double precision guesses of
`FastGaussQuadrature.gausslobatto`. Since the endpoints are known exactly, they are set
to $\mp 1$ afterwards rather than left to the root finder, so that
[`lobatto_legendre_nodes`](@ref) returns exactly `0` and `1` at the ends.

### Weights

For the Lobatto family the weights are available in closed form,

```math
w_i = \frac{2}{s \, (s-1) \, \big[ P_{s-1} (x_i) \big]^{2}} ,
```

a formula which is valid at the endpoints as well as at the interior nodes. Evaluating
$P_{s-1}$ through the same recurrence used for the nodes and rescaling to $[0,1]$
completes the rule.

```@repl rules
LobattoLegendreQuadrature(3)      # Simpson's rule
LobattoLegendreQuadrature(3)(x -> x^3)
```

As for the Gauss rule, `fast=true` selects the double precision path through
FastGaussQuadrature.jl.


## Chebyshev points

The Chebyshev-based rules all sample at points that are known in closed form, so no root
finding is needed. [`chebyshev_points`](@ref) provides the two kinds, and
[`chebyshev_nodes`](@ref) their images in $[0,1]$.

The points of the **first kind** are the roots of the Chebyshev polynomial $T_s$,

```math
x_i = \cos \left( \frac{(2i-1) \, \pi}{2s} \right) ,
\qquad i = 1, \dots, s ,
```

and lie strictly inside the interval. The implementation evaluates the algebraically
equivalent form $\sin \big( (s-2i+1) \pi / 2s \big)$, which is more accurate near the
ends of the interval, where the cosine is flat and loses relative precision.

The points of the **second kind** are the extrema of $T_{s-1}$,

```math
x_i = \cos \left( \frac{(i-1) \, \pi}{s-1} \right) ,
\qquad i = 1, \dots, s ,
```

and include both endpoints, which makes them the Chebyshev analogue of a Lobatto node
set. They require $s \ge 2$.

Both are generated in reverse index order so that the resulting vectors are ascending,
and both are evaluated in `BigFloat` arithmetic before being converted to `T`.

```@repl rules
chebyshev_points(5, 1)
chebyshev_points(5, 2)
```

### Deriving the weights

For both kinds the weights follow from expanding the integrand in a Chebyshev series
instead of the Lagrange basis. Writing $x = \cos\theta$ turns the interpolation problem
into a trigonometric one — the correspondence between Chebyshev and Fourier series that
underlies spectral methods generally [boyd2001](@cite) — the coefficients of the series
follow from the sampled values by a discrete cosine transform, and the series is integrated
term by term using

```math
\int_{-1}^{+1} T_{2j} (x) \, dx = - \frac{2}{4 j^2 - 1} ,
\qquad
\int_{-1}^{+1} T_{2j+1} (x) \, dx = 0 .
```

Only the even-order Chebyshev polynomials contribute, and each contributes a term with
denominator $4j^2-1$. This is the common origin of the weight formulae below.

### Three variants

Which node set is used distinguishes three classical rules [trefethen2008](@cite):

| nodes | rule | provided by |
|---|---|---|
| Chebyshev roots, in $(-1,1)$ | Fejér's first rule [fejer1933](@cite) | [`GaussChebyshevQuadrature`](@ref) |
| Chebyshev extrema, in $(-1,1)$ | Fejér's second rule [fejer1933](@cite) | — |
| Chebyshev extrema, in $[-1,1]$ | Clenshaw-Curtis [clenshaw1960](@cite) | [`ClenshawCurtisQuadrature`](@ref) |

The first and third are also called the *classical* and *practical* Clenshaw-Curtis
formulae. Fejér's second rule uses the extrema of $T_{s+1}$ *excluding* the endpoints and
is **not** provided by this package. Note that its weights are not interchangeable with
those of Clenshaw-Curtis: applying the Fejér-2 weight formula to the endpoint-inclusive
Chebyshev points gives a rule that is exact only for linear functions, regardless of `s`.


## Gauss-Chebyshev quadrature (Fejér's first rule)

[`GaussChebyshevQuadrature`](@ref) is the interpolatory rule on the `s` Chebyshev points
of the first kind. Carrying out the term-by-term integration described above gives

```math
w_i = \frac{2}{s} \left( 1 - 2 \sum_{j=1}^{\lfloor s/2 \rfloor}
      \frac{\cos ( 2 j \theta_i )}{4 j^2 - 1} \right) ,
      \qquad \theta_i = \frac{(2i-1) \, \pi}{2s} ,
```

with a further factor $1/2$ for the move to $[0,1]$. This rule is classically known as
**Fejér's first rule**. All its weights are positive and all its nodes are interior.

Being interpolatory on `s` nodes, it is exact for polynomials of degree $\le s-1$, so its
order is `s`. For odd `s` it gains one further degree by symmetry.

!!! warning "Not the weighted Gauss-Chebyshev rule"
    The name Gauss-Chebyshev is also used for the rule that approximates the *weighted*
    integral $\int_{-1}^{+1} f(x) \, (1-x^2)^{-1/2} \, dx$ with the equal weights
    $\pi / s$, and which is exact to degree $2s-1$. The rule implemented here shares its
    nodes with that rule but approximates the *unweighted* integral
    $\int_0^1 f(x) \, dx$, which is what a [`QuadratureRule`](@ref) evaluates. Its degree
    of exactness is correspondingly only $s-1$.

```@repl rules
quad = GaussChebyshevQuadrature(5)
order(quad)
quad(x -> x^4)
```

Its weight sum is $O(s^2)$ just as for Clenshaw-Curtis, so the same `IT=Float64` remark
applies; see [Cost](@ref) below.


## Clenshaw-Curtis quadrature

[`ClenshawCurtisQuadrature`](@ref) is the interpolatory rule on the `s` Chebyshev points
of the second kind. With $n = s-1$ and $\vartheta_k = 2 \pi k / n$, the same term-by-term
integration gives

```math
w_k = \frac{c_k}{n} \left( 1 - \sum_{j=1}^{\lfloor n/2 \rfloor}
      \frac{b_j}{4 j^2 - 1} \cos ( j \vartheta_k ) \right) ,
```

where $c_k = 1$ at the two endpoints and $2$ otherwise, and $b_j = 2$ except for the
final term of an even-length sum, where it is $1$. These two factors are the boundary
corrections of the underlying cosine transform: the endpoints are shared by only one
half-period, and the highest mode of an even-length transform is not duplicated. As
always, a factor $1/2$ maps the rule to $[0,1]$.

The rule is exact for polynomials of degree $\le s-1$, giving order `s`, with one bonus
degree for odd `s`. All weights are positive [imhof1963](@cite). It requires $s \ge 2$.

This is the explicit closed form derived in [reid2014](@citet), which is the form the
implementation follows; the rule itself goes back to [clenshaw1960](@citet). Comparing the
symbols with the code in `src/clenshaw_curtis.jl`: `n = s-1` is Reid's $N$, `c(k,n)` is
$c_k$, `b(j,n)` is $b_j$, and `ϑ(k,n)` is $\vartheta_k$.

```@repl rules
ClenshawCurtisQuadrature(3)          # Simpson's rule again
ClenshawCurtisQuadrature(9)(x -> exp(x)) - (exp(1) - 1)
```

!!! warning "The factor $b_j$ is easy to lose"
    Some presentations, including the sample tables in [reid2014](@citet), omit the factor
    $b_j = 1$ on the final term of an even-length sum — equivalently, they do not halve the
    last Chebyshev coefficient $a_N$. This matters only for even $N$, because the
    coefficients of odd order integrate to zero anyway, but there it costs a degree of
    exactness. For $N = 4$ the correct weights on $[-1,+1]$ are
    $(\tfrac{1}{15}, \tfrac{8}{15}, \tfrac{12}{15}, \tfrac{8}{15}, \tfrac{1}{15})$, exact
    to degree 5, whereas dropping the factor gives
    $(0.05, 0.5667, 0.7667, 0.5667, 0.05)$, exact only to degree 3. Both sets sum to 2, so
    the usual sanity check does not catch the difference. This package includes the factor
    and is verified against Reid's explicit formula in the test suite.

### Convergence

Because all the weights are positive and the rule is interpolatory, the error is bounded
by the best polynomial approximation error $E_n^*$ of degree $n = s-1$,

```math
| I - I_n | \le 4 \, E_n^* ,
```

which by the Weierstrass approximation theorem implies convergence for *every* continuous
integrand [trefethen2008](@cite). Gauss-Legendre satisfies the same bound with
$E_{2n+1}^*$, reflecting its doubled degree of exactness, and this is the origin of the
folklore that Clenshaw-Curtis is "half as good".

In practice it is not — an observation reported as early as [ohara1968](@citet), though it
took a long time to become widely known. For an integrand whose $k$-th derivative has
bounded variation, Clenshaw-Curtis obeys the *same* algebraic bound as Gauss, with $2n$
rather than $n$ [trefethen2008](@cite):

```math
| I - I_n | \le \frac{32 \, V}{15 \, \pi \, k \, (2n+1-k)^{k}} .
```

The mechanism is aliasing. On the Chebyshev grid $T_{n+p}$ and $T_{n-p}$ are
indistinguishable, so the rule returns $I(T_{n-p})$ when handed $T_{n+p}$; since
$I(T_{n-p})$ is itself $O(n^{-2})$ small, the error contributed by the first Chebyshev
coefficients beyond the exactness limit is far smaller than a naive count of exact degrees
suggests. Gauss quadrature has a decisive advantage only when $f$ is analytic in a sizable
neighbourhood of the interval, where it converges like $\rho^{-2n}$ against
$\rho^{-n}$ — and there both methods reach machine precision so quickly that the
difference rarely matters.

The following reproduces the comparison of [trefethen2008](@citet), Figure 2, for three
integrands of decreasing smoothness. The last column is the ratio of the Clenshaw-Curtis
error to the Gauss-Legendre error at the same number of nodes:

```@example rules
using Printf
symmetric(quad, f) = 2 * quad(ξ -> f(2ξ - 1))     # [0,1] rule applied on [-1,+1]

#                                                  exact value of ∫₋₁¹ f(x) dx
cases = [("1/(1+16x^2)", x -> 1/(1+16x^2),         atan(4)/2),
         ("exp(-1/x^2)", x -> x == 0 ? zero(x) : exp(-1/x^2),
                                                   2*(exp(-1) - sqrt(π)*0.15729920705028513)),
         ("|x|^3",       x -> abs(x)^3,            0.5)]

for (name, f, exact) in cases, n in (8, 16, 24)
    eg = abs(symmetric(GaussLegendreQuadrature(n+1; fast=true), f) - exact)
    ec = abs(symmetric(ClenshawCurtisQuadrature(n+1; IT=Float64), f) - exact)
    @printf("%-12s n=%-3d Gauss %.2e   Clenshaw-Curtis %.2e   ratio %5.2f\n",
            name, n, eg, ec, ec/eg)
end
```

For $1/(1+16x^2)$, analytic but with poles at $\pm i/4$ close to the interval, the ratio
sits at a little over two. For the two non-analytic integrands it hovers around one, and
for $\exp(-1/x^2)$ Clenshaw-Curtis is at times the more accurate of the two. Individual
entries are noisy — the $n = 24$ row for $\exp(-1/x^2)$ catches the Gauss error at a
particularly favourable point — so it is the trend rather than any single ratio that
matters. Nowhere does the asymptotic factor of two in the degree of exactness translate
into a factor of two in accuracy.

That factor becomes visible only for polynomials and entire functions: $x^{20}$ is
integrated exactly by Gauss from $n \ge 10$ but by Clenshaw-Curtis only from $n \ge 20$,
and for $e^x$ both reach machine precision well before $n = 16$.

### Cost

The weight sum above is evaluated directly, at a cost of $O(s^2)$ operations, and by
default in `BigFloat`. Passing `IT=Float64` is about two orders of magnitude faster and is
the appropriate choice whenever the result is wanted in `Float64`:

```@repl rules
ClenshawCurtisQuadrature(Float64, 64; IT=Float64) ≈ ClenshawCurtisQuadrature(64)
```

Computing the weights through a fast cosine transform instead reduces this to
$O(s \log s)$; see [gentleman1972](@citet) and [waldvogel2006](@citet). That is not done
here, because the package's priority is arbitrary-precision accuracy for moderate node
counts rather than throughput at large $s$.


## Lobatto-Chebyshev quadrature

[`LobattoChebyshevQuadrature`](@ref) is the interpolatory rule on the Chebyshev points of
the second kind — but those are exactly the Clenshaw-Curtis nodes, and an interpolatory
rule is uniquely determined by its nodes. The two rules are therefore identical, and the
implementation delegates accordingly:

```@repl rules
LobattoChebyshevQuadrature(6) == ClenshawCurtisQuadrature(6)
lobatto_chebyshev_points(6) == clenshaw_curtis_points(6)
```

The separate name is retained because it is the natural counterpart to
[`GaussChebyshevQuadrature`](@ref) within the Chebyshev family, mirroring the
Gauss/Lobatto distinction of the Legendre family.


## Umbrella constructor

[`ChebyshevQuadrature`](@ref) selects between the two Chebyshev rules by the kind of the
underlying points, in the same way that [`chebyshev_points`](@ref) does:

```@repl rules
ChebyshevQuadrature(4, 1) == GaussChebyshevQuadrature(4)
ChebyshevQuadrature(4, 2) == LobattoChebyshevQuadrature(4)
```


## Summary

| rule | nodes | endpoints included | order | requires |
|---|---|---|---|---|
| [`RiemannQuadratureLeft`](@ref) | $0$ | one | 1 | |
| [`RiemannQuadratureRight`](@ref) | $1$ | one | 1 | |
| [`MidpointQuadrature`](@ref) | $1/2$ | no | 2 | |
| [`TrapezoidalQuadrature`](@ref) | $0, 1$ | both | 2 | |
| [`GaussLegendreQuadrature`](@ref) | roots of $P_s$ | no | $2s$ | |
| [`LobattoLegendreQuadrature`](@ref) | roots of $P_{s-1}'$ and $\pm 1$ | both | $2s-2$ | $s \ge 2$ |
| [`GaussChebyshevQuadrature`](@ref) | Chebyshev, first kind | no | $s$ | |
| [`ClenshawCurtisQuadrature`](@ref) | Chebyshev, second kind | both | $s$ | $s \ge 2$ |
| [`LobattoChebyshevQuadrature`](@ref) | Chebyshev, second kind | both | $s$ | $s \ge 2$ |

The orders quoted for the Chebyshev-based rules are guaranteed values; for an odd number
of nodes these rules integrate one additional degree exactly.

Not provided: Radau rules, which fix one endpoint, and Fejér's second rule, the
interpolatory rule on the Chebyshev extrema *excluding* the endpoints.
