```@meta
CurrentModule = QuadratureRules
```

# Numerical Quadrature

The need for numerical quadrature often arises for evaluating the definite integral of a function that has no explicit antiderivative (indefinite integral) or whose antiderivative is not easy to obtain.
The idea of numerical quadrature is to approximate $\int_{a}^{b} f(x) \, dx$ by a sum $\sum_{i=1}^{s} b_{i} f(x_{i})$.

The points $x_i$ at which the integrand is sampled are called the *nodes* of the rule and
the coefficients $b_i$ its *weights*. A quadrature rule is thus completely described by
these two vectors, together with a statement of how accurate the resulting approximation is.
This is precisely the information stored in a [`QuadratureRule`](@ref).


## The reference interval

Every rule in this package is defined on the reference interval $[0,1]$, so that

```math
\int_0^1 f(x) \, dx \approx \sum_{i=1}^{s} b_i \, f(c_i) ,
\qquad c_i \in [0,1] .
```

An integral over an arbitrary interval $[a,b]$ is recovered by the affine change of
variables $x = a + (b-a) \, \xi$,

```math
\int_a^b f(x) \, dx = (b-a) \int_0^1 f \big( a + (b-a) \, \xi \big) \, d\xi
\approx (b-a) \sum_{i=1}^{s} b_i \, f \big( a + (b-a) \, c_i \big) .
```

Because $\int_0^1 1 \, dx = 1$, the weights of every rule on the reference interval sum
to one. This normalisation is what makes rules from different families directly
comparable, and it is checked throughout the test suite.

Most of the classical theory is instead formulated on the symmetric interval $[-1,+1]$,
which is where orthogonal polynomials such as the Legendre and Chebyshev polynomials are
defined. The package follows this convention internally and offers both views:

- the `*_points` and `*_point_weights` functions work on $[-1,+1]$,
- the `*_nodes` and `*_weights` functions work on $[0,1]$,

related by $c_i = (x_i + 1) / 2$, with the weights scaled correspondingly by a factor
$1/2$ so that they sum to $1$ rather than $2$. See
[Points, Nodes and Weights](@ref) for the full list.


## Interpolatory quadrature

Almost all useful quadrature rules are *interpolatory*: given $s$ distinct nodes $c_i$,
one replaces the integrand by the unique polynomial $p$ of degree $s-1$ that interpolates
$f$ at those nodes, and integrates that polynomial exactly. Writing the interpolant in
the Lagrange basis,

```math
\ell_i (x) = \prod_{j \neq i} \frac{x - c_j}{c_i - c_j} ,
\qquad
p(x) = \sum_{i=1}^{s} f(c_i) \, \ell_i (x) ,
```

and integrating term by term identifies the weights as the integrals of the basis
functions,

```math
b_i = \int_0^1 \ell_i (x) \, dx .
```

Two consequences are worth stating explicitly, because they explain much of what follows.
First, once the nodes are fixed, the weights are determined; there is exactly one
interpolatory rule per node set. This is why [`LobattoChebyshevQuadrature`](@ref) and
[`ClenshawCurtisQuadrature`](@ref), which share their nodes, are necessarily the same
rule. Second, an interpolatory rule on $s$ nodes reproduces every polynomial of degree
$\le s-1$ exactly, since such a polynomial is its own interpolant.


## Degree of exactness and order

A rule has *degree of exactness* $d$ if it integrates every polynomial of degree $\le d$
exactly but fails for some polynomial of degree $d+1$. Since integration is linear, it is
enough to test the monomials, which on $[0,1]$ gives the *moment conditions*

```math
\sum_{i=1}^{s} b_i \, c_i^{k} = \int_0^1 x^k \, dx = \frac{1}{k+1} ,
\qquad k = 0, 1, \dots, d .
```

Throughout this package the *order* $p$ reported by [`order`](@ref) is related to the
degree of exactness by $p = d + 1$: a rule of order $p$ integrates polynomials of degree
$\le p - 1$ exactly. The order is the quantity that governs the convergence rate of a
composite rule, and it is the convention used for the order of the collocation and
variational integrators that this package is primarily written for.

The order stored in a rule is *sharp*: the rule really does fail for some polynomial of
degree $p$. In particular, the Clenshaw-Curtis and Chebyshev rules pick up one extra
degree of exactness when the number of nodes is odd — the additional monomial they would
fail on is odd about the midpoint of the interval and therefore integrates to zero on
both sides — and the reported order includes that bonus.

Sharpness has a useful consequence. Since the order is a function of the nodes and
weights rather than of the family a rule was constructed from, two rules that coincide
compare equal. The three-node Clenshaw-Curtis and Lobatto-Legendre rules are both
Simpson's rule, and both report order 4.


## How accurate can a rule be?

An $s$-node rule has $2s$ free parameters, namely the $s$ nodes and the $s$ weights.
Matching $2s$ moment conditions therefore suggests a maximal degree of exactness of
$2s-1$, and this bound is attained. The rules that attain it are the *Gauss* rules, and
they are characterised by a remarkable property: the nodes must be the roots of the
degree-$s$ polynomial that is orthogonal to all polynomials of lower degree with respect
to the inner product

```math
\langle f, g \rangle = \int_{-1}^{+1} f(x) \, g(x) \, w(x) \, dx .
```

For the unweighted integral, $w \equiv 1$, these are the Legendre polynomials $P_s$, and
the resulting rule is [`GaussLegendreQuadrature`](@ref). To see why orthogonality is the
right condition, write an arbitrary polynomial $f$ of degree $\le 2s-1$ as
$f = q \, P_s + r$ with $q$ and $r$ of degree $\le s-1$. The rule integrates $r$ exactly
because it is interpolatory, and the remaining term $\int q \, P_s$ vanishes by
orthogonality while its quadrature sum vanishes because $P_s$ is zero at every node.

Two further properties follow from this construction and matter in practice: all Gauss
nodes lie in the open interval, and all Gauss weights are positive. Positive weights mean
that the rule does not amplify rounding errors in the sampled function values, which is
why high-order Gauss rules remain numerically well behaved where interpolation through
equidistant points does not.


## Constrained node families

Sometimes it is useful to prescribe some of the nodes in advance, at the cost of
accuracy. Each node fixed in advance removes one free parameter and hence one degree of
exactness:

| family | prescribed nodes | free nodes | degree of exactness | order |
|---|---|---|---|---|
| Gauss | none | $s$ | $2s-1$ | $2s$ |
| Radau | one endpoint | $s-1$ | $2s-2$ | $2s-1$ |
| Lobatto | both endpoints | $s-2$ | $2s-3$ | $2s-2$ |

The Lobatto family is implemented as [`LobattoLegendreQuadrature`](@ref). Including the
endpoints is worth two degrees of exactness whenever the values of the integrand at the
endpoints are needed anyway, as in finite element methods, collocation schemes and
variational integrators, where the endpoint values are shared between neighbouring
elements or time steps.

The Radau family is implemented as [`RadauLegendreQuadrature`](@ref). It prescribes one
endpoint and leaves the other free, which is the asymmetric middle ground between the two:
a single point of contact with the boundary, bought for a single degree of exactness. Which
endpoint is prescribed is a genuine choice rather than a convention, so it is passed
explicitly as `:left` or `:right`. The one-sided constraint is what makes these nodes the
basis of the stiffly accurate Radau IA and Radau IIA collocation methods.


## Chebyshev-based rules

Gauss rules require the roots of an orthogonal polynomial, which have no closed form and
must be computed numerically. An alternative is to fix the nodes at points that are known
analytically and accept the lower degree of exactness $s-1$ that any interpolatory rule
on $s$ nodes provides. The natural choice is the Chebyshev points, that is, the
projections onto the interval of equally spaced points on the unit semicircle:

```math
x_i = \cos \left( \frac{(2i-1) \, \pi}{2s} \right)
\quad \text{(first kind)} ,
\qquad
x_i = \cos \left( \frac{(i-1) \, \pi}{s-1} \right)
\quad \text{(second kind)} .
```

Points of the first kind lie strictly inside the interval; points of the second kind
include both endpoints, making them the natural "Lobatto-like" choice. Both cluster
towards the ends of the interval, which is exactly what makes polynomial interpolation
through them stable, in contrast to interpolation through equidistant points.

The weights are obtained by expanding the integrand in a Chebyshev series rather than in
the Lagrange basis. The coefficients of that series follow from the sampled values by a
discrete cosine transform, and the series can be integrated term by term using

```math
\int_{-1}^{+1} T_{2j}(x) \, dx = - \frac{2}{4 j^2 - 1} ,
\qquad
\int_{-1}^{+1} T_{2j+1}(x) \, dx = 0 ,
```

which is the origin of the $4j^2 - 1$ denominators appearing in the weight formulae of
[`GaussChebyshevQuadrature`](@ref) and [`ClenshawCurtisQuadrature`](@ref).

Although their degree of exactness is only about half that of a Gauss rule with the same
number of nodes, these rules converge at a comparable rate for smooth integrands. The
reason is that convergence for non-polynomial integrands is governed by how well the
function is approximated by polynomials on the node set, not by degree of exactness
alone, and Chebyshev interpolation is near-optimal in that respect.

Two results make this precise [trefethen2008](@cite). First, because the Chebyshev weights
are positive [imhof1963](@cite) and the rules are interpolatory, the quadrature error is
bounded by $4 E_n^*$, where $E_n^*$ is the error of the best polynomial approximation of
degree $n$; a Gauss rule satisfies the same bound with $E_{2n+1}^*$. Second, and less
obviously, for integrands of finite smoothness Clenshaw-Curtis satisfies the same
*algebraic* error bound as Gauss, with $2n$ in place of $n$.

The mechanism behind the second result is aliasing. On a grid of Chebyshev points the
polynomials $T_{n+p}$ and $T_{n-p}$ take identical values, so the rule cannot tell them
apart and returns $I(T_{n-p})$ when handed $T_{n+p}$. Since those exact integrals are
themselves $O(n^{-2})$ small, the Chebyshev coefficients just past the exactness limit
contribute far less error than merely counting exact degrees would suggest. A Gauss rule,
by contrast, is exact up to degree $2n+1$ and then fails abruptly. The upshot is that
Gauss quadrature has a decisive advantage only for integrands analytic in a sizable
neighbourhood of the interval — and there both methods converge geometrically, so the
distinction is of little practical consequence. See the
[Convergence](@ref) discussion for a numerical comparison.


## Variable transformations and the double-exponential formula

Every rule discussed so far is interpolatory: it fixes a set of nodes, integrates the
polynomial through them, and is characterised by a degree of exactness. There is a second,
quite different way to build an accurate rule, and [`TanhSinhQuadrature`](@ref) is the one
representative of it in this package.

The starting point is that the humble trapezoidal rule, which is only of order 2 on a finite
interval, becomes *spectrally* accurate on the whole real line. If $g$ is analytic and
suitably small in the strip $|\operatorname{Im} t| < d$, the infinite trapezoidal sum with
step $h$ satisfies

```math
\Bigg| \int_{-\infty}^{\infty} g(t) \, dt - h \sum_{k=-\infty}^{\infty} g(k h) \Bigg|
= \mathcal{O} \big( e^{-2 \pi d / h} \big) ,
```

because the error is controlled by the Fourier transform of $g$, which decays
exponentially at the rate set by the width of the strip [trefethen2014](@cite). The
difficulty is that the sum is infinite. It can only be truncated cheaply if $g$ itself
decays very fast.

That is exactly what a well-chosen change of variables arranges. Substituting

```math
x = \tanh \left( \frac{\pi}{2} \sinh t \right)
```

maps $\mathbb{R}$ onto the open interval $(-1,+1)$, and

```math
\frac{dx}{dt} = \frac{\tfrac{\pi}{2} \cosh t}{\cosh^2 \big( \tfrac{\pi}{2} \sinh t \big)}
= \mathcal{O} \big( e^{- \tfrac{\pi}{2} e^{|t|}} \big) .
```

The transformed integrand $f(x(t)) \, x'(t)$ therefore vanishes *double
exponentially* — hence the name **double-exponential formula**, due to
[takahasi1974](@citet) — and the sum may be truncated after a few dozen terms. Balancing the
discretisation error against the truncation error gives, for $N$ nodes,

```math
| I - I_N | = \mathcal{O} \big( e^{-c N / \log N} \big) ,
```

which is optimal in a precise sense among transformations of this kind [mori2001](@cite);
[mori2005](@citet) recounts how the formula was found.

Two features distinguish such a rule sharply from the interpolatory ones.

First, **it has no degree of exactness at all.** The truncated sum reproduces no polynomial
exactly, not even a constant: the weights sum to one only up to the truncation error. The
order reported by [`order`](@ref) is consequently `0` for tanh-sinh, which is the honest
statement that no polynomial is integrated exactly, and its accuracy has to be described by
the convergence rate above instead. In practice each halving of $h$ roughly *doubles* the
number of correct digits.

Second, the transformation sends both endpoints to infinity, so the nodes never reach them
and an integrand that is singular there can be passed in unchanged. This is what tanh-sinh
is for: no interpolatory rule copes with $\int_0^1 x^{-1/2} \, dx$, because polynomial
approximation of $x^{-1/2}$ near the origin is hopeless, whereas tanh-sinh handles it as a
matter of course.

The catch is that the benefit is bounded by the arithmetic. A node of type `T` cannot
approach an endpoint closer than about `eps(T)`, so for an integrand growing like
$x^{-1/2}$ the neglected tail is of size $\sqrt{\texttt{eps(T)}}$ — about 8 digits in
`Float64` — however many levels are used. Increasing the *precision* rather than the number
of nodes is what buys accuracy here, which is why tanh-sinh is the standard tool of
high-precision numerical integration [bailey2005](@cite) and sits naturally in a package
built around arbitrary-precision rules.


## Composite rules and convergence

For a fixed rule of order $p$, subdividing $[a,b]$ into $n$ subintervals of length
$h = (b-a)/n$ and applying the rule on each gives a composite rule with error

```math
\Bigg| \int_a^b f(x) \, dx
     - h \sum_{j=1}^{n} \sum_{i=1}^{s} b_i \, f \big( a + (j-1) \, h + h \, c_i \big) \Bigg|
= \mathcal{O} \big( h^{p} \big)
```

for sufficiently smooth $f$. This is the sense in which order is the practically relevant
number: doubling the number of subintervals reduces the error by a factor $2^{p}$.

The alternative strategy is to keep a single interval and increase the number of nodes.
For analytic integrands this converges geometrically rather than algebraically, which is
why a Gauss or Clenshaw-Curtis rule with a few dozen nodes routinely reaches machine
precision on a smooth integrand.


## Arbitrary precision

The nodes of the Gauss and Lobatto rules are roots of polynomials and are only available
numerically. Computing them in double precision limits the resulting rule to roughly
`Float64` accuracy, which is not enough when the rule is used to construct a
high-order integrator whose coefficients must satisfy order conditions to full precision.

This package therefore computes nodes and weights in an internal working precision,
controlled by the keyword argument `IT` and defaulting to `BigFloat`, and converts the
result to the requested element type `T` only at the very end. The roots themselves are
obtained by taking the double precision approximations from
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) as
initial guesses and refining them with Newton's method in the working precision. The
closed-form node and weight formulae of the Chebyshev family are likewise evaluated in `IT`,
so that all intermediate quantities, and not merely the final result, are computed in high
precision.

This matters more than it might appear. Setting `IT` equal to `T` — computing a `Float64`
rule entirely in `Float64` — is much faster, but the round-off accumulated in the
intermediate terms leaves the moment conditions violated by a few hundred times `eps(T)`,
whereas computing in `BigFloat` and rounding only at the end satisfies them to better than
`eps(T)`, i.e. returns correctly rounded nodes and weights. Since the whole point of the
package is to supply coefficients that meet the order conditions of a high-order integrator
to full precision, `BigFloat` is the default everywhere.

Where the extra accuracy is not needed, the Legendre rules accept `fast=true`, which
takes the nodes and weights directly from FastGaussQuadrature.jl in double precision.


## Exact and symbolic arithmetic

Everything said so far describes what happens for a numeric element type. High precision is a
workaround for arithmetic that rounds, though, and if `T` does not round there is nothing to
work around: for an element type from outside the numeric tower the working precision defaults
to `T` itself, and nodes and weights come out exactly.

This is what happens with the element type of a computer algebra system such as
[SymPy.jl](https://github.com/JuliaPy/SymPy.jl) or
[SymPyPythonCall.jl](https://github.com/jverzani/SymPyPythonCall.jl):

```julia
julia> using QuadratureRules, SymPyPythonCall

julia> simplify.(gauss_legendre_nodes(typeof(Sym(1)), 2))
2-element Vector{Sym{PythonCall.Py}}:
 1/2 - sqrt(3)/6
 sqrt(3)/6 + 1/2

julia> simplify.(radau_legendre_weights(typeof(Sym(1)), 3, :right))
3-element Vector{Sym{PythonCall.Py}}:
 4/9 - sqrt(6)/36
 sqrt(6)/36 + 4/9
             1/9
```

Two things change on this path, and neither is visible from the outside. The roots of the
Legendre polynomial can no longer be refined from a double precision guess, so they are
computed exactly instead, as the eigenvalues of the polynomial's companion matrix — which a
computer algebra system resolves into radicals. And the Chebyshev and Clenshaw-Curtis closed
forms are evaluated at the exact ``\pi`` rather than at a rounded one. The closed forms for
the weights are already written in terms of the element type and need no change at all.

Downstream this is what lets a package like
[RungeKutta.jl](https://github.com/JuliaGNI/RungeKutta.jl) tabulate its methods symbolically,
by asking for the coefficients at a symbolic element type and typesetting the result.

Three caveats. Exact node computation is only as good as the algebra system behind it; SymPy
resolves the roots as radicals up to the quartic and falls back on implicit representations
beyond it, so this is a facility for the small stage numbers that tableaus use, not a
replacement for the numerical path. It also needs `eigvals` for the element type, which is
what a computer algebra system supplies and what the exact members of the numeric tower do
not: `Rational` and `Complex` element types are therefore counted as numeric and computed in
`BigFloat` like the floating point ones, since a quadrature node is an algebraic number that
no rational can represent to begin with. And [`TanhSinhQuadrature`](@ref) does not take part:
its node count is not given in advance but follows from where nodes and weights stop being
resolvable in `T`, which a type that does not round cannot answer, so it accepts floating
point element types only and throws an `ArgumentError` for anything else.

For an element type that merely *wraps* numbers rather than computing exactly, pass `IT`
explicitly — `gauss_legendre_nodes(T, s; IT=BigFloat)` computes the rule in `BigFloat` as
usual and converts at the end, exactly as it does for a numeric `T`.

```@example
using QuadratureRules

quad = GaussLegendreQuadrature(BigFloat, 5)
sum(weights(quad) .* nodes(quad).^9) - 1/big(10)   # exact to full precision
```
