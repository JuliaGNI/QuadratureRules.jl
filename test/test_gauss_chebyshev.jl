import FastTransforms
import QuadratureRules: shift_nodes, unshift_nodes

@testset "$(rpad("Gauß-Chebyshev",80))" begin

    for s in 2:10
        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(s, 1)
        @test GaussChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(1))
        
        @test sum(weights(GaussChebyshevQuadrature(s))) ≈ 1

        x = reverse(FastTransforms.chebyshevpoints(Float64, s, Val(1)))
        c = shift_nodes(x)

        @test chebyshev_points(Float64, s, Val(1)) ≈ x
        @test chebyshev_nodes(Float64, s, Val(1))  ≈ c
        @test nodes(GaussChebyshevQuadrature(s))   ≈ c

        @test gauss_chebyshev_points(Float64, s) == chebyshev_points(Float64, s, Val(1))
        @test gauss_chebyshev_nodes(Float64, s)  == chebyshev_nodes(Float64, s, Val(1))
        @test gauss_chebyshev_points(s) == gauss_chebyshev_points(Float64, s)
        @test gauss_chebyshev_nodes(s)  == gauss_chebyshev_nodes(Float64, s)
        @test chebyshev_points(s, 1) == chebyshev_points(Float64, s, Val(1))
        @test chebyshev_nodes(s, 1)  == chebyshev_nodes(Float64, s, Val(1))

        @test unshift_nodes(chebyshev_nodes(Float64, s, Val(1))) ≈ chebyshev_points(Float64, s, Val(1))

        # Gauss-Chebyshev is Fejér's first rule: an interpolatory rule on the s
        # Chebyshev points of the first kind, hence exact up to degree s-1
        @test order(GaussChebyshevQuadrature(s)) == s

        let q = GaussChebyshevQuadrature(BigFloat, s)
            @test all(0 ≤ cᵢ ≤ 1 for cᵢ in nodes(q))
            @test issorted(nodes(q), lt = <)
            @test all(bᵢ > 0 for bᵢ in weights(q))

            for k in 0:order(q)-1
                @test sum(weights(q) .* nodes(q).^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
            end
        end

    end

end
