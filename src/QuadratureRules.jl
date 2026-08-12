module QuadratureRules

    import FastGaussQuadrature
    import Polynomials
    import Polynomials: Polynomial

    include("utils.jl")


    include("quadrature_rule.jl")

    export QuadratureRule
    export nnodes, nodes, order, weights

    include("tabulated_quadratures.jl")

    export RiemannQuadratureLeft,
           RiemannQuadratureRight,
           MidpointQuadrature,
           TrapezoidalQuadrature

    include("chebyshev.jl")
    include("clenshaw_curtis.jl")
    include("gauss_legendre.jl")
    include("lobatto_legendre.jl")

    export ChebyshevQuadrature,
           GaussChebyshevQuadrature,
           LobattoChebyshevQuadrature,
           ClenshawCurtisQuadrature,
           GaussLegendreQuadrature,
           LobattoLegendreQuadrature

    export chebyshev_points,
           gauss_chebyshev_points,
           lobatto_chebyshev_points,
           clenshaw_curtis_points,
           gauss_legendre_points,
           lobatto_legendre_points

    export chebyshev_nodes,
           gauss_chebyshev_nodes,
           lobatto_chebyshev_nodes,
           clenshaw_curtis_nodes,
           gauss_legendre_nodes,
           lobatto_legendre_nodes


end
