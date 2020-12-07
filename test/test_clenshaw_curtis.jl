import FastTransforms
import QuadratureRules: scale_weights, shift_nodes

@testset "$(rpad("Clenshaw-Curtis",80))" begin

    @test_throws ErrorException ClenshawCurtisQuadrature(1)

    for s in 2:10
        @test ClenshawCurtisQuadrature(s) == ClenshawCurtisQuadrature(Float64, s)

        @test sum(weights(GaussChebyshevQuadrature(s))) ≈ 1
        
        μ = FastTransforms.chebyshevmoments1(Float64, s)
        b = scale_weights(reverse(FastTransforms.clenshawcurtisweights(μ)))
        c = shift_nodes(reverse(FastTransforms.clenshawcurtisnodes(Float64, s)))

        @test nodes(ClenshawCurtisQuadrature(s)) ≈ c
        @test weights(ClenshawCurtisQuadrature(s)) ≈ b
    end

end
