using QuadratureRules
using Test

include("test_utils.jl")
include("test_quadrature_rule.jl")
include("test_gauss_chebyshev.jl")
include("test_gauss_legendre.jl")
include("test_lobatto_chebyshev.jl")
include("test_lobatto_legendre.jl")
include("test_radau_legendre.jl")
include("test_clenshaw_curtis.jl")
include("test_tanh_sinh.jl")
include("test_tabulated_quadratures.jl")
include("test_order.jl")
include("test_precision.jl")

# test_symbolic.jl is deliberately not included here. It needs SymPyPythonCall, which brings a
# private Python installation with it, so it is not part of this package's test target but of
# test/symbolic/Project.toml, and runs on its own:
#
#     julia --project=test/symbolic -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
#     julia --project=test/symbolic test/test_symbolic.jl
