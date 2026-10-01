using ExplicitImports
using QuadratureRules
using Test

# Implicit imports, stale explicit imports, explicit imports and qualified accesses of names
# that are not public or through a module that does not own the name, and self-qualified
# accesses fail this test.
test_explicit_imports(QuadratureRules)
