# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

### K1 · The comment of `test/symbolic/Project.toml` compares it with `docs/Project.toml`, which resolves the package differently

- location: `test/symbolic/Project.toml:9`
- evidence: the comment says that `Pkg.develop(PackageSpec(path=pwd()))` adds QuadratureRules
  "from the checkout, as it does for docs/Project.toml". `docs/Project.toml` has
  `[sources] QuadratureRules = {path = ".."}`, so the docs environment does not need that
  `Pkg.develop`.
- kind: docs
- found: 2026-10-01
