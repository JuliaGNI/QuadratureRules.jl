```@meta
CurrentModule = QuadratureRules
```

# References

The theory summarised in [Numerical Quadrature](@ref) and the derivations in
[Quadrature Rules](@ref) follow the sources listed below.

For the Chebyshev-based rules, [trefethen2008](@citet) is the natural starting point: it
settles the question of how Clenshaw-Curtis compares with Gauss quadrature in practice and
gives the taxonomy of the three Chebyshev variants. The construction of the Clenshaw-Curtis
weights used in this package follows [reid2014](@citet); [waldvogel2006](@citet) gives the
fast $O(n \log n)$ alternative. [davis1984](@citet) is the standard monograph on numerical
integration, and [trefethen2019](@citet) covers the approximation theory that the convergence
results rest on.

```@bibliography
```
