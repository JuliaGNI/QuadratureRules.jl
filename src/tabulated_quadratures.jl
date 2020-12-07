
"""
Riemannian quadrature on the left.
"""
function RiemannQuadratureLeft(T=Float64)
    QuadratureRule(1, T[0], T[1])
end

"""
Riemannian quadrature on the right.
"""
function RiemannQuadratureRight(T=Float64)
    QuadratureRule(1, T[1], T[1])
end

"""
Midpoint quadrature.
"""
function MidpointQuadrature(T=Float64)
    QuadratureRule(2, T[1//2], T[1])
end

"""
Trapezoidal quadrature.
"""
function TrapezoidalQuadrature(T=Float64)
    QuadratureRule(2, T[0, 1], T[1//2, 1//2])
end
