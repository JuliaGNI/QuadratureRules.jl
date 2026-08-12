module QuadratureRules

    import FastGaussQuadrature
    import Polynomials
    import Polynomials: Polynomial

    # extended below rather than defined here, so that the packages of the ecosystem share
    # one generic function per accessor instead of one per package
    import GeometricBase: nnodes, nodes, order, weights

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
    include("radau_legendre.jl")
    include("tanh_sinh.jl")

    export ChebyshevQuadrature,
           GaussChebyshevQuadrature,
           LobattoChebyshevQuadrature,
           ClenshawCurtisQuadrature,
           GaussLegendreQuadrature,
           LobattoLegendreQuadrature,
           RadauLegendreQuadrature,
           TanhSinhQuadrature

    export symmetric_chebyshev_nodes,
           symmetric_gauss_chebyshev_nodes,
           symmetric_lobatto_chebyshev_nodes,
           symmetric_clenshaw_curtis_nodes,
           symmetric_gauss_legendre_nodes,
           symmetric_lobatto_legendre_nodes,
           symmetric_radau_legendre_nodes,
           symmetric_tanh_sinh_nodes

    export chebyshev_nodes,
           gauss_chebyshev_nodes,
           lobatto_chebyshev_nodes,
           clenshaw_curtis_nodes,
           gauss_legendre_nodes,
           lobatto_legendre_nodes,
           radau_legendre_nodes,
           tanh_sinh_nodes

    export symmetric_chebyshev_weights,
           symmetric_gauss_chebyshev_weights,
           symmetric_lobatto_chebyshev_weights,
           symmetric_clenshaw_curtis_weights,
           symmetric_gauss_legendre_weights,
           symmetric_lobatto_legendre_weights,
           symmetric_radau_legendre_weights,
           symmetric_tanh_sinh_weights

    export chebyshev_weights,
           gauss_chebyshev_weights,
           lobatto_chebyshev_weights,
           clenshaw_curtis_weights,
           gauss_legendre_weights,
           lobatto_legendre_weights,
           radau_legendre_weights,
           tanh_sinh_weights

end
