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
an arbitrary number of nodes `s` and in arbitrary precision. Both are interpolatory.
[`TanhSinhQuadrature`](@ref) stands outside this division and is treated last: it is not
interpolatory but the trapezoidal rule after a change of variables, and it is parameterised
by a refinement level rather than by a number of nodes.


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
Instead of differentiating the Legendre polynomial, the implementation forms

```math
D(x) = \frac{d^{\,s-2}}{dx^{\,s-2}} \, \big( 1 - x^2 \big)^{s-1}
```

directly with `Polynomials.derivative`. By Rodrigues' formula for the Jacobi polynomials,

```math
P^{(1,1)}_{s-2} (x) \; \propto \; \big( 1 - x^2 \big)^{-1} \,
    \frac{d^{\,s-2}}{dx^{\,s-2}} \, \big( 1 - x^2 \big)^{s-1} ,
```

and $P_{s-1}' \propto P^{(1,1)}_{s-2}$, so $D$ has the $s-2$ interior Lobatto points among
its roots — *and*, from the factor $(1-x^2)$ it retains, the two endpoints as well. $D$ has
degree $s$ and its $s$ roots are therefore exactly the $s$ Lobatto points, which is why the
whole node set is recovered from this single polynomial.

The roots are again Newton-refined from the double precision guesses of
`FastGaussQuadrature.gausslobatto`, all `s` of them at once. Since the endpoints are known
exactly, they are set
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


## Radau-Legendre quadrature

[`RadauLegendreQuadrature`](@ref) constrains exactly one endpoint of the interval to be a
node, so it sits between the Gauss rules, which constrain none, and the Lobatto rules,
which constrain both. With $s-1$ free nodes the degree of exactness is $2s-2$ and the order
$2s-1$. Unlike the Lobatto rules it is defined for `s == 1`, where the single node is the
prescribed endpoint carrying the whole weight, so the rule degenerates into a Riemann sum.

Which endpoint is prescribed changes the rule, so it is not defaulted but passed explicitly
as `:left` (the node $-1$, or $0$ after rescaling) or `:right` (the node $+1$, or $1$). The
two variants are mirror images of one another. They are the node families underlying the
Radau IA and Radau IIA collocation methods respectively; the asymmetry is the point, since
prescribing the right endpoint is what makes an implicit Runge-Kutta method stiffly
accurate. The rules go back to [radau1880](@citet); see [gautschi2000](@citet) for the
modern treatment of the Jacobi-weighted case and [hairer1996](@citet) for their use in the
numerical solution of stiff and differential-algebraic equations.

### Nodes

The left Radau points are the $s$ roots of

```math
R(x) = P_{s-1} (x) + P_s (x) ,
```

one of which is exactly $-1$, since $P_k(-1) = (-1)^k$ makes the two terms cancel there. To
see that the remaining $s-1$ roots are the right ones, note that

```math
\frac{P_{s-1} (x) + P_s (x)}{1 + x} \; \propto \; P^{(0,1)}_{s-1} (x) ,
```

the Jacobi polynomial for the weight $1-x$ on $[-1,+1]$. Its roots are precisely the free
nodes of the Gauss rule for that weight, which is what remains after the factor $1+x$ has
absorbed the prescribed endpoint — exactly the construction
`FastGaussQuadrature.gaussradau` performs, and the reason $R$ recovers the whole node set
from a single polynomial, as $D$ does in the Lobatto case.

The roots are Newton-refined from the double precision guesses of
`FastGaussQuadrature.gaussradau` and the prescribed endpoint is then set to $-1$ exactly, so
that [`radau_legendre_nodes`](@ref) returns exactly `0` there. The right points are obtained
by reflection, $x \mapsto -x$, which makes the two variants exact mirror images rather than
two independent root finds.

```@repl rules
radau_legendre_points(3, :left)
radau_legendre_points(3, :right)
```

### Weights

The weights, too, are available in closed form,

```math
w_i = \frac{1 \mp x_i}{s^2 \, \big[ P_{s-1} (x_i) \big]^{2}} ,
```

with the upper sign for the left variant and the lower one for the right. As in the Lobatto
case the single formula covers the prescribed endpoint as well: there
$\big[ P_{s-1} (\mp 1) \big]^2 = 1$, so it collapses to the familiar $2/s^2$.

```@repl rules
RadauLegendreQuadrature(2, :right)
RadauLegendreQuadrature(3, :right)(x -> x^4)      # exact, degree 4 = 2s-2
RadauLegendreQuadrature(3, :right)(x -> x^5)      # not exact, degree 5
RadauLegendreQuadrature(1, :left) == RiemannQuadratureLeft()
```

As for the Gauss and Lobatto rules, `fast=true` selects the double precision path through
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
and both are evaluated in the working precision `IT`, `BigFloat` by default, before being
converted to `T`.

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

Being interpolatory on `s` nodes, it is exact for polynomials of degree $\le s-1$. For
odd `s` it gains one further degree, the monomial of degree `s` being odd about the
midpoint of the interval, so its order is `s` for even `s` and `s+1` for odd `s`. With one
node it reduces to the midpoint rule.

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

The rule is exact for polynomials of degree $\le s-1$, with one bonus degree for odd `s`,
so its order is `s` for even `s` and `s+1` for odd `s`. All weights are positive
[imhof1963](@cite). It requires $s \ge 2$. With three nodes it is Simpson's rule and
therefore equal, order included, to the three-node Lobatto-Legendre rule.

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
default in `BigFloat`. Lowering the working precision with `IT=Float64` is one to two
orders of magnitude faster:

```@repl rules
ClenshawCurtisQuadrature(Float64, 64; IT=Float64) ≈ ClenshawCurtisQuadrature(64)
```

`BigFloat` is nevertheless the default, and deliberately so. The weights are a sum of
$O(s)$ cosine terms, and every term contributes round-off; evaluating the closed forms in
`BigFloat` and rounding only at the end delivers nodes and weights correct to the full
precision of `T`, whereas a lower working precision merely comes close. This is what the
package is for — the coefficients of a high-order integrator have to satisfy their order
conditions to full precision, and a rule that is a few units in the last place off will not
do. Lower `IT` only when the cost matters and that accuracy does not.

Computing the weights through a fast cosine transform instead reduces the cost to
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


## Tanh-Sinh quadrature

[`TanhSinhQuadrature`](@ref) is the odd one out. It is not an interpolatory rule on a
prescribed node set but the trapezoidal rule applied after a change of variables, the
*double-exponential formula* of [takahasi1974](@citet). See
[Variable transformations and the double-exponential formula](@ref) for why that is a good
idea at all; this section covers how it is computed.

### The substitution

Setting

```math
x = \tanh \left( \frac{\pi}{2} \sinh t \right) , \qquad t \in \mathbb{R} ,
```

turns $\int_{-1}^{+1} f(x) \, dx$ into an integral over the whole real line, to which the
trapezoidal rule with step $h$ is applied at the abscissae $t = k h$:

```math
\int_{-1}^{+1} f(x) \, dx \approx h \sum_{k \in \mathbb{Z}} w_k \, f(x_k) ,
\qquad
x_k = \tanh \left( \frac{\pi}{2} \sinh (k h) \right) ,
\qquad
w_k = \frac{\tfrac{\pi}{2} \cosh (k h)}{\cosh^2 \big( \tfrac{\pi}{2} \sinh (k h) \big)} .
```

The weights are just $dx/dt$ evaluated at the abscissae. Since $dx/dt$ decays like
$\exp ( - \tfrac{\pi}{2} e^{|t|} )$, the sum can be truncated after a modest number of
terms.

### The logistic form

Mapping to $[0,1]$ turns the transformation into the logistic sigmoid,

```math
c_k = \frac{1 + x_k}{2} = \frac{1}{1 + e^{-\pi \sinh (k h)}} ,
\qquad
1 - c_k = \frac{1}{1 + e^{+\pi \sinh (k h)}} ,
```

and this, not $(1 + \tanh u)/2$, is what the implementation evaluates. The reason is
cancellation: the whole point of the outermost nodes is *how close they are to the
endpoints*, and $1 - \tanh u$ loses all its significant digits for large $u$, whereas
$1/(1 + e^{2u})$ delivers the small number directly. The rule is assembled from these
offsets in symmetric pairs about the centre node $c_0 = 1/2$, so that the weight vector
comes out symmetric to the last bit and the nodes ascending by construction.

### Levels and nesting

The parameter of the rule is not a node count but a **level** `n`, which fixes the step size
$h = 2^{-n}$. Halving $h$ retains every previous abscissa and inserts one new one between
each pair, so the levels are nested — a property inherited from the trapezoidal rule and the
reason the classical implementations refine level by level, reusing all previous integrand
values:

```@repl rules
issubset(tanh_sinh_nodes(2), tanh_sinh_nodes(3))
nnodes.(TanhSinhQuadrature.(1:5))
```

Because the truncation criterion below is a threshold in $t$ that all levels share, the node
set of level `n` is simply the multiples of $2^{-n}$ below that threshold, and the
implementation can generate it in a single loop at the finest step instead.

The nesting holds while the rule is still resolving, which is the range worth using anyway.
Past it the pairwise merge described below folds the outermost nodes of consecutive levels
together, and the inclusion fails — from level 5 in `Float32` and from level 6 in `Float64`.

### Truncation

Where to stop is decided by the target type `T`, not by a tolerance. The sum is truncated at
the first `k` whose weight rounds to zero in `T` or whose node rounds to an endpoint — of
either $[0,1]$ or $[-1,+1]$, since the latter carries one bit less next to the endpoints and
would otherwise degenerate. Should a pair still be indistinguishable in `T` from its
predecessor, its weight is folded into that predecessor rather than added as a new node,
which leaves the quadrature sum untouched and keeps the nodes strictly increasing.

The upshot is that **no node ever coincides with an endpoint**, which is precisely what
allows an integrand singular there to be handed over as it stands:

```@repl rules
quad = TanhSinhQuadrature(3);
extrema(nodes(quad))
quad(x -> 1 / sqrt(x * (1 - x))) - π      # ∫₀¹ dx/√(x(1-x)) = π
```

### No degree of exactness

The truncated sum integrates no polynomial exactly, not even a constant — its weights sum to
one only up to the truncation error — so [`order`](@ref) reports `0`. Tanh-sinh is the only
rule in this package for which the order says nothing about the accuracy:

```@repl rules
order(TanhSinhQuadrature(3))
sum(weights(TanhSinhQuadrature(BigFloat, 2))) - 1
```

### Doubling of digits per level

What describes the accuracy instead is the convergence rate
$\mathcal{O} ( e^{-cN/\log N} )$, which in practice means that each level roughly doubles the
number of correct digits until the precision of `T` is exhausted:

```@example rules
using Printf

#                                                    exact value of ∫₀¹ f(x) dx
cases = [("exp(x)",     exp,                          exp(big(1)) - 1),
         ("log(x)",     log,                          big(-1)),
         ("√x·log(1/x)", x -> sqrt(x) * log(1/x),     big(4)//9)]

for (name, f, exact) in cases
    errs = [Float64(abs(TanhSinhQuadrature(BigFloat, n)(f) - exact)) for n in 1:5]
    @printf("%-12s %s\n", name, join((@sprintf("%9.2e", e) for e in errs), " "))
end
```

The first column is level 1 and the last level 5; the final entries are at the `BigFloat`
round-off level and no longer measure the rule. Note that the logarithmic singularity at
$x = 0$ costs nothing at all — it is integrated just as accurately as $e^x$.

### Comparison with Gauss-Legendre

At an equal number of nodes the two rules are good at opposite things:

```@example rules
for (level, N) in ((1, 13), (3, 51))
    ts = TanhSinhQuadrature(level)
    gl = GaussLegendreQuadrature(N; fast=true)

    for (name, f, exact) in [("exp(x)",    exp,            exp(1) - 1),
                             ("log(x)",    log,            -1.0),
                             ("1/sqrt(x)", x -> 1/sqrt(x),  2.0)]
        @printf("%2d nodes  %-12s tanh-sinh %8.2e   Gauss-Legendre %8.2e\n",
                N, name, abs(ts(f) - exact), abs(gl(f) - exact))
    end
end
```

For an integrand analytic up to and including the endpoints, Gauss-Legendre is far ahead: at
13 nodes it has already reached machine precision on $e^x$, where tanh-sinh is still at
$10^{-5}$, because 13 nodes buy Gauss exactness to degree 25 while tanh-sinh has spent most of
its nodes resolving the ends of an interval where nothing interesting happens. As soon as
there is an endpoint singularity the comparison reverses, and not by a small margin.
Gauss-Legendre is left with about four correct digits on $\log x$ and barely two on
$x^{-1/2}$ even at 51 nodes, and adding nodes helps it only algebraically, because polynomial
approximation of those integrands near the origin is hopeless. Tanh-sinh reaches machine
precision on the first and its intrinsic limit on the second.

The rule of thumb follows: use [`GaussLegendreQuadrature`](@ref) unless the integrand
misbehaves at an endpoint, and tanh-sinh when it does.

### The `√eps(T)` limit

That limit deserves a closer look, because it is the one thing that must be understood before
relying on the rule. A node of type `T` cannot come closer to an endpoint than about
`eps(T)`, so the tail of the transformed integrand beyond the last node is not negligible if
$f$ grows there. For $f(x) = x^{-1/2}$ it is of size $\sqrt{\texttt{eps(T)}}$, and no
further level helps:

```@example rules
for T in (Float32, Float64, BigFloat)
    errs = [Float64(abs(TanhSinhQuadrature(T, n)(x -> 1/sqrt(x)) - 2)) for n in 3:6]
    @printf("%-9s %s   √eps = %8.2e\n", T,
            join((@sprintf("%8.2e", e) for e in errs), " "), sqrt(eps(T)))
end
```

The remedy is not more nodes but more precision — 8 correct digits in `Float64`, 39 at the
default `BigFloat` precision, 78 at twice that. This is the mechanism behind the use of
tanh-sinh as the workhorse of high-precision quadrature [bailey2005](@cite), and the
clearest illustration of what this package computes in arbitrary precision *for*.

A logarithmic singularity, by contrast, is integrated to full precision at any type, as the
convergence table above shows: the tail it leaves behind is of size
$\texttt{eps(T)} \, \log \texttt{eps(T)}$, which is negligible.

### The useful range of levels

Each level doubles the number of nodes, so the node count grows like $2^n$ while the
attainable accuracy stops improving once the precision of `T` is reached. There is no point
in going beyond level 3 for `Float64` or level 5 at the default `BigFloat` precision, and the
tables above are the evidence. Constructing the rule costs one `sinh`, `cosh` and `exp` per
node in the working precision `IT`; unlike the Chebyshev rules there is no summation, so
`IT=BigFloat` is comparatively cheap and there is little reason to lower it.


## Summary

| rule | nodes | endpoints included | order | requires |
|---|---|---|---|---|
| [`RiemannQuadratureLeft`](@ref) | $0$ | one | 1 | |
| [`RiemannQuadratureRight`](@ref) | $1$ | one | 1 | |
| [`MidpointQuadrature`](@ref) | $1/2$ | no | 2 | |
| [`TrapezoidalQuadrature`](@ref) | $0, 1$ | both | 2 | |
| [`GaussLegendreQuadrature`](@ref) | roots of $P_s$ | no | $2s$ | |
| [`LobattoLegendreQuadrature`](@ref) | roots of $P_{s-1}'$ and $\pm 1$ | both | $2s-2$ | $s \ge 2$ |
| [`RadauLegendreQuadrature`](@ref) | roots of $P_{s-1} + P_s$ | one | $2s-1$ | |
| [`GaussChebyshevQuadrature`](@ref) | Chebyshev, first kind | no | $s$ / $s+1$ | |
| [`ClenshawCurtisQuadrature`](@ref) | Chebyshev, second kind | both | $s$ / $s+1$ | $s \ge 2$ |
| [`LobattoChebyshevQuadrature`](@ref) | Chebyshev, second kind | both | $s$ / $s+1$ | $s \ge 2$ |
| [`TanhSinhQuadrature`](@ref) | $\tanh \big( \tfrac{\pi}{2} \sinh (k h) \big)$, $h = 2^{-n}$ | no | — | $n \ge 1$ |

The two orders quoted for the Chebyshev-based rules are for an even and an odd number of
nodes respectively; these rules integrate one additional degree exactly when `s` is odd.
Every order in the table is sharp, so rules that coincide agree in their order too:
`ClenshawCurtisQuadrature(3) == LobattoLegendreQuadrature(3)`,
`GaussChebyshevQuadrature(1) == MidpointQuadrature()` and
`RadauLegendreQuadrature(1, :left) == RiemannQuadratureLeft()`.

The nodes of the Radau rule are stated for the left variant; the right variant is its
reflection, and which one is meant is selected by the `endpoint` argument rather than
defaulted.

The dash in the last row is not an omission. Tanh-sinh is the one rule here with no degree of
exactness whatever, so [`order`](@ref) reports `0` and its accuracy is described by the
convergence rate discussed above instead. It also differs in taking a level `n` rather than a
node count `s`; the number of nodes follows from the truncation and grows like $2^n$.

Not provided: Fejér's second rule, the interpolatory rule on the Chebyshev extrema
*excluding* the endpoints.
