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

# The symbolic tests need SymPyPythonCall, which brings a private Python installation with it.
# It is therefore not part of the package's test target but of test/symbolic/Project.toml, and
# the tests run whenever the active environment provides it, i.e. in that environment.
if isnothing(Base.identify_package("SymPyPythonCall"))
    @info "SymPyPythonCall is not available in this environment, skipping the symbolic tests. " *
          "Run them with --project=test/symbolic."
else
    include("test_symbolic.jl")
end
