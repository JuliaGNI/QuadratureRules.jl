using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Utility functions" include("utils.jl")
    @safetestset "Quadrature rule" include("quadrature_rule.jl")
    @safetestset "Gauss-Chebyshev" include("gauss_chebyshev.jl")
    @safetestset "Gauss-Legendre" include("gauss_legendre.jl")
    @safetestset "Lobatto-Chebyshev" include("lobatto_chebyshev.jl")
    @safetestset "Lobatto-Legendre" include("lobatto_legendre.jl")
    @safetestset "Radau-Legendre" include("radau_legendre.jl")
    @safetestset "Clenshaw-Curtis" include("clenshaw_curtis.jl")
    @safetestset "Tanh-sinh" include("tanh_sinh.jl")
    @safetestset "Tabulated quadrature rules" include("tabulated_quadratures.jl")
    @safetestset "Order" include("integration/order.jl")
    @safetestset "Working precision" include("integration/precision.jl")
end
if "slow" in GROUPS
    @safetestset "Doctests" include("quality/doctests.jl")
    @safetestset "Symbolic" include("test_symbolic.jl")
end
