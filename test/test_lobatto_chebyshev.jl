import FastTransforms
import QuadratureRules: shift_nodes

@testset "$(rpad("Lobatto-Chebyshev",80))" begin

    @test_throws ErrorException LobattoChebyshevQuadrature(1)

    for s in 2:10
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(s, 2)
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(2))
        
        @test sum(weights(LobattoChebyshevQuadrature(s))) ≈ 1

        c = shift_nodes(reverse(FastTransforms.chebyshevpoints(Float64, s, Val(2))))

        @test chebyshev_points(Float64, s, Val(2)) ≈ c
        @test nodes(LobattoChebyshevQuadrature(s)) ≈ c

        # TODO test weights

    end

end
