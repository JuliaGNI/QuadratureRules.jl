# Changelog

All notable changes to QuadratureRules.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the usual
Julia reading of `0.x`: a bump of the *minor* version may break.

## [Unreleased]

### Added

- **The nodes and weights are now documented on the unit interval as well as the symmetric
  one.** Every closed form was stated only in the `SymmetricInterval` variable ``x_i``, even
  though `UnitInterval` is the default, so a reader who knows a node family from the
  Runge-Kutta literature — where it is given on ``[0,1]`` — could not match it against the
  text. The `*_nodes` docstrings and the manual now also give the shifted characterisations
  ``P_s(2x-1)`` (Gauss), ``d^{s-2}/dx^{s-2}((x-x^2)^{s-1})`` (Lobatto) and
  ``d^{s-1}/dx^{s-1}(x^s(x-1)^{s-1})`` / ``d^{s-1}/dx^{s-1}(x^{s-1}(x-1)^s)`` (Radau left and
  right), and the `*_weights` docstrings the corresponding closed forms in the unit-interval
  nodes ``c_i``. The Lobatto form is in fact simpler there, the halving cancelling the
  numerator. The Radau expressions are documentation only: the implementation still reflects
  the left nodes rather than evaluating the mirrored polynomial, so that the two variants stay
  exact mirror images.
- **The Rodrigues formula for the Legendre polynomial** is documented alongside the three-term
  recurrence, in both `QuadratureRules._legendre` and the manual, with the reason the
  recurrence is used for evaluation instead. It is what identifies the derivative form of the
  Lobatto nodes as an antiderivative of ``P_{s-1}``, a step the manual previously asserted
  without stating the identity it rests on. The Rodrigues formula for the Jacobi polynomial
  ``P^{(0,1)}_{s-1}`` is now stated in the Radau section for the same reason: it is what the
  unit-interval Radau polynomials rest on, and it makes explicit that the differentiated
  product is its *numerator*, the division by ``1+x`` being what removes the prescribed
  endpoint.
- **The moment conditions are described as a construction, not only as a test.** For a
  prescribed node set they are a transposed-Vandermonde system in the weights; the manual now
  says so, notes that this package uses closed forms instead because that system is badly
  conditioned for larger `s`, and mentions that the Runge-Kutta literature defines several
  families through exactly these conditions under the name of the simplifying assumption
  ``B(s)``.

### Fixed

- **`QuadratureRules._legendre` now carries the three-term recurrence upwards in a loop.** It
  recursed on ``j-1`` and ``j-2``, so it re-evaluated the same ``P_k`` exponentially often —
  a degree ``30`` polynomial cost almost three million calls — which made building the node
  polynomial the dominant cost for larger `s`. The arithmetic is unchanged and performed in
  the same order, so every node and weight the package returns is identical bit for bit.

## [0.2.0] – 2026-08-13

The first release since the package was substantially reworked. It is **breaking**: the
accessors are renamed, three rules report a different order, one returns different weights,
and four functions are now GeometricBase generics. Each change is listed below with the
reason and, where numbers moved, what they moved from and to.

### ⚠️ Breaking

- **The interval is now a keyword argument, not part of the function name.** The `*_points`
  and `*_point_weights` functions are gone. Each family provides one node function and one
  weight function, both taking an `interval` keyword argument:

  ```julia
  gauss_legendre_nodes(2)                                    # [0,1], as before
  gauss_legendre_nodes(2; interval = SymmetricInterval())     # was gauss_legendre_points(2)
  gauss_legendre_weights(2; interval = SymmetricInterval())   # was gauss_legendre_point_weights(2)
  ```

  `interval` takes the new exported singletons `UnitInterval()` for ``[0,1]``, the default and
  the interval of `QuadratureRule`, or `SymmetricInterval()` for ``[-1,+1]``, both subtypes of
  the new abstract type `QuadratureInterval`. What used to be encoded by two different words —
  *point* versus *node*, neither of which names an interval — is now an argument that says
  which interval it is. `SymmetricInterval` describes the interval, not the node set: a Radau
  rule has deliberately asymmetric nodes on either interval. The keyword is annotated
  `::QuadratureInterval`, so anything else — a bare `:symmetric`, say — is a `TypeError`
  naming the keyword, never a silent fallback to `[0,1]`.

  This drops 16 exported names and matches how the package already selects the Radau endpoint
  and the Chebyshev kind. The 16 removed functions map to keyword calls mechanically:
  `symmetric_X_nodes(args...)` becomes `X_nodes(args...; interval = SymmetricInterval())`, and
  likewise for the weights. There is no deprecation path; the old names were public for a
  single patch release, and `*_point_weights` was never released at all.

- **`LobattoChebyshevQuadrature` returns different weights.** It combined two different
  point sets: its nodes were the Chebyshev points of the second kind, which include the
  endpoints, but its weights were the Fejér-2 weights for the *interior* points
  ``\cos(i\pi/(s+1))``. The resulting rule was exact only for linear functions, at every
  `s`. For `s = 3` it returned the weights `[1/3, 1/3, 1/3]` on the nodes `[0, 1/2, 1]`,
  where the correct interpolatory weights are Simpson's `[1/6, 2/3, 1/6]`. It now delegates
  to `ClenshawCurtisQuadrature`, with which it shares its nodes.

- **The reported `order` changed for three rules.** `GaussChebyshevQuadrature` implements
  Fejér's first rule for the unweighted integral; its weights were correct, but it reported
  order `2s`, which is the exactness of the *weighted* Gauss-Chebyshev rule for
  ``\int f(x)/\sqrt{1-x^2} \, dx`` rather than of the rule actually evaluated. The sharp
  orders are now reported:

  | rule | before | after |
  | --- | --- | --- |
  | `GaussChebyshevQuadrature` | `2s` | `isodd(s) ? s+1 : s` |
  | `LobattoChebyshevQuadrature` | `2s` | `isodd(s) ? s+1 : s` |
  | `ClenshawCurtisQuadrature` | `s` | `isodd(s) ? s+1 : s` |

  The Legendre orders are unchanged (`2s`, `2s-2`, and `2s-1` for the new Radau rules).
  `TanhSinhQuadrature` reports `order == 0` by design, being exact for no polynomial.

- **`nnodes`, `nodes`, `order` and `weights` are now GeometricBase generics.** They are
  imported from `GeometricBase` and extended rather than defined here, so that the packages
  of the ecosystem share one generic function per accessor instead of one per package. They
  are still exported. Code that defined its own methods on `QuadratureRules.nodes` now
  extends a different function object. This adds a hard dependency on
  `GeometricBase = "0.14.8"`.

- **The default working precision now depends on the element type.** The `IT` keyword
  defaults to `_default_arithmetic(T)` instead of a hard-coded `BigFloat`: `BigFloat` for
  the numeric tower (`AbstractFloat`, `Rational`, `Integer`, `Complex`) and `T` itself for
  anything else, on the assumption that such a type does its own arithmetic exactly.
  Behaviour for `Float64` and `BigFloat` is unchanged; element types outside the numeric
  tower now compute in themselves. `IT=BigFloat` recovers the old behaviour for them.
  `TanhSinhQuadrature` keeps a literal `IT=BigFloat`, since it needs a type that rounds.

- **The `@big` macro and the `_big` helpers are removed** from `src/utils.jl`, replaced by
  threading the working precision `IT` explicitly. They were unexported, but were listed
  among the documented internals.

### Added

- `QuadratureInterval`, `UnitInterval` and `SymmetricInterval`, described above.

- **Tanh-sinh quadrature**: `TanhSinhQuadrature` together with `tanh_sinh_nodes` and
  `tanh_sinh_weights`. The rule for integrands with endpoint singularities, converging
  double-exponentially. It differs from the other families in three ways: it takes a *level*
  `n` rather than a node count, the number of nodes following from where they cease to be
  resolvable in `T`; it is computed on ``[0,1]`` in the logistic form, which yields the
  outermost nodes without cancellation; and it accepts floating point element types only,
  throwing an `ArgumentError` for anything else.

- **Radau-Legendre quadrature**: `RadauLegendreQuadrature` together with
  `radau_legendre_nodes` and `radau_legendre_weights`. Prescribes one endpoint
  and leaves the rest free, for order `2s-1`. Which endpoint is prescribed is a further
  positional argument, `:left` (Radau IA) or `:right` (Radau IIA), deliberately without a
  default, as the two variants are not interchangeable. The right variant is obtained by
  reflecting the left one, which makes the two exact mirror images of each other.

- **Weights accessors for every family**, i.e. the eight `*_weights` functions. Previously
  only the nodes were available without constructing the rule.

- **Nodes and weights on symbolic element types.** For an element type outside the numeric
  tower the roots are computed exactly, as the eigenvalues of the companion matrix, and the
  closed forms are evaluated in the type itself, ``\pi`` included. A symbolic `T` therefore
  yields nodes and weights in closed form — `gauss_legendre_nodes(typeof(Sym(1)), 2)` gives
  `1/2 ± sqrt(3)/6` — rather than `BigFloat` values wrapped in symbolic clothing. Exercised
  by a separate CI job against SymPyPythonCall, in the environment of
  `test/symbolic/Project.toml`.

- An `IT` keyword throughout, including `GaussChebyshevQuadrature`, which previously ignored
  the keywords it was given.

- `unscale_weights`, the inverse of `scale_weights` (unexported).

- **Documentation**: new pages on the theory of quadrature (`quadrature.md`), on the rules
  one by one (`rules.md`) and a bibliography (`references.md`, via `DocumenterCitations`).
  The Clenshaw-Curtis derivation is documented with its sources and its convergence, and
  verified against the literature.

- New test suites for the working precision across `(T, IT)` combinations, for the reported
  orders, and for the Radau, tanh-sinh and symbolic paths.

### Changed

- The Legendre nodes are computed once per rule construction rather than twice, the
  constructors sharing the closed form for the weights with the accessors.

### Fixed

- `ClenshawCurtisQuadrature` honours its `IT` keyword. The `@big` macro inside the weight
  sum forced `BigFloat` regardless of `IT`, so the weights came out `BigFloat` while the
  nodes followed `IT`, and the constructor rejected the mismatched element types: every
  `IT` other than `BigFloat` threw a `MethodError`. `LobattoChebyshevQuadrature`, which
  forwards `IT`, was affected in the same way.
- Several docstring and manual corrections: the sign of the node pair in the tanh-sinh back
  end, the Jacobi weight of the Radau nodes, the Lobatto weights cross-references, and two
  over-broad claims about tanh-sinh exactness.

## [0.1.10] – 2026-08-12

### Added

- `*_points` and `*_nodes` accessors for all rules, so that nodes can be obtained without
  constructing a rule. (Folded into the `interval` keyword in the unreleased changes above.)

### Changed

- The `fast` keyword is accepted by the quadrature constructors only, no longer by the node
  functions.

## [0.1.9] – 2026-08-12

### Changed

- Simplified the quadrature evaluation and the computation of the Clenshaw-Curtis weights.
- Julia 1.10 is now the minimum supported version; compat entries updated.

## [0.1.8] – 2026-08-12

### Removed

- The `GenericLinearAlgebra` dependency.

### Changed

- Added `.JuliaFormatter.toml`; updated the GitHub Actions workflows to `setup-julia` v3.

## [0.1.7] – 2026-03-27

### Changed

- CI runs on Julia 1.11, 1.12 and 1.13; support for 1.6 dropped.
- `GenericLinearAlgebra 0.4` compat.

## [0.1.6] – 2023-10-09

### Fixed

- A typo and a bug in the Clenshaw-Curtis quadrature.

### Changed

- `FastGaussQuadrature 1` compat.

## [0.1.5] – 2023-09-20

### Changed

- `Polynomials 4` compat; explicit Julia compat entry; CI and Documenter updates.

## [0.1.4] – 2022-11-18

### Changed

- `FastGaussQuadrature 0.5` compat; added the Register workflow.

## [0.1.3] – 2022-03-09

### Changed

- Julia 1.4 is now the minimum supported version.
- `Polynomials 2`, then `3`, and `GenericLinearAlgebra 0.3` compat; CI on Julia 1.6–1.8.

## [0.1.2] – 2020-12-17

### Removed

- The `SpecialPolynomials` dependency.

## [0.1.1] – 2020-12-15

### Added

- The `GenericLinearAlgebra` dependency.

## [0.1.0] – 2020-12-07

Initial release: the `QuadratureRule` type and its functor, the tabulated rules (left and
right Riemann, midpoint, trapezoidal), and the Clenshaw-Curtis, Chebyshev and Legendre
rules.

[Unreleased]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.10...v0.2.0
[0.1.10]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.9...v0.1.10
[0.1.9]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.8...v0.1.9
[0.1.8]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.7...v0.1.8
[0.1.7]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.6...v0.1.7
[0.1.6]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.4...v0.1.5
[0.1.4]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/JuliaGNI/QuadratureRules.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JuliaGNI/QuadratureRules.jl/releases/tag/v0.1.0
