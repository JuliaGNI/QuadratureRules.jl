using Documenter
using QuadratureRules

# the same doctest setup as docs/make.jl
DocMeta.setdocmeta!(QuadratureRules, :DocTestSetup, :(using QuadratureRules); recursive = true)

# the manual pages resolve `CurrentModule = QuadratureRules` in `Main`, and a
# `@safetestset` runs this file in a module of its own
@eval Main import QuadratureRules

doctest(QuadratureRules)
