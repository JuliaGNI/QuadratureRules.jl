module QuadratureRules

import FastGaussQuadrature
import GenericLinearAlgebra
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
include("tanh_sinh.jl")

export ChebyshevQuadrature,
    GaussChebyshevQuadrature,
    LobattoChebyshevQuadrature,
    ClenshawCurtisQuadrature,
    GaussLegendreQuadrature,
    LobattoLegendreQuadrature,
    TanhSinhQuadrature

end
