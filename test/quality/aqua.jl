using Aqua
using QuadratureRules
using Test

# Package-level quality assurance: type piracy, method ambiguities, stale and duplicated
# dependencies, undefined exports, unbound type parameters, `Project.toml` validity.
Aqua.test_all(QuadratureRules)
