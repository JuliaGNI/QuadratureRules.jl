import FastTransforms
import QuadratureRules: chebyshev_points, shift_nodes

@testset "$(rpad("Gauß-Chebyshev",80))" begin

    for s in 2:10
        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(s, 1)
        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(1))
        
        @test sum(weights(GaussChebyshevQuadrature(s))) ≈ 1

        c = shift_nodes(reverse(FastTransforms.chebyshevpoints(Float64, s, Val(1))))

        @test chebyshev_points(Float64, s, Val(1)) ≈ c
        @test nodes(GaussChebyshevQuadrature(s))   ≈ c

        # TODO test weights

    end

end
