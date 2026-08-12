@testset "$(rpad("Tabulated Quadrature Rules",80))" begin

    @test typeof(RiemannQuadratureLeft()) <: QuadratureRule
    @test typeof(RiemannQuadratureRight()) <: QuadratureRule
    @test typeof(MidpointQuadrature()) <: QuadratureRule
    @test typeof(TrapezoidalQuadrature()) <: QuadratureRule

    for quad in (RiemannQuadratureLeft, RiemannQuadratureRight,
                 MidpointQuadrature, TrapezoidalQuadrature)
        @test sum(weights(quad())) == 1
        @test eltype(quad()) == Float64
        @test eltype(quad(BigFloat)) == BigFloat
        @test nnodes(quad()) == length(nodes(quad()))
        @test all(0 ≤ cᵢ ≤ 1 for cᵢ in nodes(quad()))

        # a rule of order p integrates polynomials up to degree p-1 exactly
        let q = quad(BigFloat)
            for k in 0:order(q)-1
                @test sum(weights(q) .* nodes(q).^k) == 1 / BigFloat(k+1)
            end
        end
    end

    @test nodes(RiemannQuadratureLeft())    == [0.0]
    @test weights(RiemannQuadratureLeft())  == [1.0]
    @test order(RiemannQuadratureLeft())    == 1

    @test nodes(RiemannQuadratureRight())   == [1.0]
    @test weights(RiemannQuadratureRight()) == [1.0]
    @test order(RiemannQuadratureRight())   == 1

    @test nodes(MidpointQuadrature())       == [0.5]
    @test weights(MidpointQuadrature())     == [1.0]
    @test order(MidpointQuadrature())       == 2

    @test nodes(TrapezoidalQuadrature())    == [0.0, 1.0]
    @test weights(TrapezoidalQuadrature())  == [0.5, 0.5]
    @test order(TrapezoidalQuadrature())    == 2

end
