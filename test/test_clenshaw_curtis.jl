import FastTransforms
import QuadratureRules: scale_weights, shift_nodes, unshift_nodes

@testset "$(rpad("Clenshaw-Curtis",80))" begin

    @test_throws ErrorException ClenshawCurtisQuadrature(1)

    for s in 2:10
        @test ClenshawCurtisQuadrature(s) == ClenshawCurtisQuadrature(Float64, s)

        @test sum(weights(ClenshawCurtisQuadrature(s))) ≈ 1
        @test order(ClenshawCurtisQuadrature(s)) == s

        # an s-node interpolatory rule integrates polynomials up to degree s-1 exactly
        let q = ClenshawCurtisQuadrature(BigFloat, s)
            for k in 0:order(q)-1
                @test sum(weights(q) .* nodes(q).^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
            end
        end


        μ = FastTransforms.chebyshevmoments1(Float64, s)
        b = scale_weights(reverse(FastTransforms.clenshawcurtisweights(μ)))
        c = shift_nodes(reverse(FastTransforms.clenshawcurtisnodes(Float64, s)))

        @test nodes(ClenshawCurtisQuadrature(s)) ≈ c
        @test weights(ClenshawCurtisQuadrature(s)) ≈ b

        @test clenshaw_curtis_nodes(Float64, s) ≈ c

        # the quadrature computes its nodes in IT = BigFloat and rounds to T,
        # whereas clenshaw_curtis_nodes(Float64, s) computes them in Float64,
        # so the two agree exactly only at equal working precision
        @test clenshaw_curtis_nodes(Float64, s)  ≈  nodes(ClenshawCurtisQuadrature(s))
        @test clenshaw_curtis_nodes(BigFloat, s) == nodes(ClenshawCurtisQuadrature(BigFloat, s))
        @test clenshaw_curtis_points(Float64, s) ≈ reverse(FastTransforms.clenshawcurtisnodes(Float64, s))
        @test unshift_nodes(clenshaw_curtis_nodes(Float64, s)) ≈ clenshaw_curtis_points(Float64, s)

        @test clenshaw_curtis_points(s) == clenshaw_curtis_points(Float64, s)
        @test clenshaw_curtis_nodes(s)  == clenshaw_curtis_nodes(Float64, s)
    end

end
