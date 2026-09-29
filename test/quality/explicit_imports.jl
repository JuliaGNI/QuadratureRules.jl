using ExplicitImports
using QuadratureRules
using Test

# Stale explicit imports, explicit imports and qualified accesses through a module that does
# not own the name, and self-qualified accesses fail this test.
test_explicit_imports(
    QuadratureRules;
    # QuadratureRules has no implicit import; the check is off because this guard
    # keeps the one form that the ecosystem packages share, and that form does not require it
    no_implicit_imports = false,
    # `nnodes`, `nodes`, `order` and `weights` are imported from GeometricBase, which does
    # not declare them `public`; that is an open GeometricBase decision
    all_explicit_imports_are_public = false,
    # every qualified access is public; the check is off for the same shared form, so
    # that a dependency which drops a `public` declaration does not fail this package's tests
    all_qualified_accesses_are_public = false
)
