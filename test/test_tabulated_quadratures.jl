@testset "$(rpad("Tabulated Quadrature Rules",80))" begin

    @test typeof(RiemannQuadratureLeft()) <: QuadratureRule
    @test typeof(RiemannQuadratureRight()) <: QuadratureRule
    @test typeof(MidpointQuadrature()) <: QuadratureRule
    @test typeof(TrapezoidalQuadrature()) <: QuadratureRule

end
