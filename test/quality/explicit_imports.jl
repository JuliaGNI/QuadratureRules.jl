using ExplicitImports
using QuadratureRules
using Test

# Stale explicit imports, qualified accesses to names their module does not own, and
# self-qualified accesses fail this test.
test_explicit_imports(
    QuadratureRules;
    # QuadratureRules has no implicit import today; the check is off because this guard
    # keeps the one form that the ecosystem packages share, and that form does not require it
    no_implicit_imports = false,
    # `nnodes`, `nodes`, `order` and `weights` are imported from GeometricBase, which does
    # not declare them `public`; that is an open GeometricBase decision
    all_explicit_imports_are_public = false,
    # every qualified access is public today; the check is off for the same shared form, so
    # that a dependency which drops a `public` declaration does not fail this package's tests
    all_qualified_accesses_are_public = false
)
