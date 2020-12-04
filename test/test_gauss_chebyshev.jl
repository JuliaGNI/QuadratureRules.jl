import FastTransforms
import QuadratureRules: shift_nodes

@testset "$(rpad("Gauß-Chebyshev",80))" begin

    for s in 2:10
        c = shift_nodes(reverse(FastTransforms.chebyshevpoints(Float64, s, Val(1)), 1))

        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(s, Val(1))
        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(1))
        
        @test nodes(GaussChebyshevQuadrature(s)) ≈ c

        # TODO test weights

    end

end
