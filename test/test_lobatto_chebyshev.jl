import FastTransforms
import QuadratureRules: shift_nodes, unshift_nodes

@testset "$(rpad("Lobatto-Chebyshev",80))" begin

    @test_throws ErrorException LobattoChebyshevQuadrature(1)
    @test_throws ErrorException chebyshev_points(Float64, 1, Val(2))
    @test_throws ErrorException chebyshev_nodes(Float64, 1, Val(2))
    @test_throws ErrorException lobatto_chebyshev_points(1)
    @test_throws ErrorException lobatto_chebyshev_nodes(1)

    for s in 2:10
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(s, 2)
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(2))
        
        @test sum(weights(LobattoChebyshevQuadrature(s))) ≈ 1

        x = reverse(FastTransforms.chebyshevpoints(Float64, s, Val(2)))
        c = shift_nodes(x)

        @test chebyshev_points(Float64, s, Val(2)) ≈ x
        @test chebyshev_nodes(Float64, s, Val(2))  ≈ c
        @test nodes(LobattoChebyshevQuadrature(s)) ≈ c

        @test lobatto_chebyshev_points(Float64, s) == chebyshev_points(Float64, s, Val(2))
        @test lobatto_chebyshev_nodes(Float64, s)  == chebyshev_nodes(Float64, s, Val(2))
        @test lobatto_chebyshev_points(s) == lobatto_chebyshev_points(Float64, s)
        @test lobatto_chebyshev_nodes(s)  == lobatto_chebyshev_nodes(Float64, s)
        @test chebyshev_points(s, 2) == chebyshev_points(Float64, s, Val(2))
        @test chebyshev_nodes(s, 2)  == chebyshev_nodes(Float64, s, Val(2))

        @test unshift_nodes(chebyshev_nodes(Float64, s, Val(2))) ≈ chebyshev_points(Float64, s, Val(2))

        # the Chebyshev points of the second kind are exactly the Clenshaw-Curtis
        # nodes, so the interpolatory rule on them is the Clenshaw-Curtis rule
        @test LobattoChebyshevQuadrature(s) == ClenshawCurtisQuadrature(s)
        @test LobattoChebyshevQuadrature(BigFloat, s) == ClenshawCurtisQuadrature(BigFloat, s)
        @test order(LobattoChebyshevQuadrature(s)) == (isodd(s) ? s+1 : s)

        let q = LobattoChebyshevQuadrature(BigFloat, s)
            @test nodes(q)[begin] == 0
            @test nodes(q)[end]   == 1
            @test issorted(nodes(q), lt = <)
            @test all(bᵢ > 0 for bᵢ in weights(q))

            for k in 0:order(q)-1
                @test sum(weights(q) .* nodes(q).^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
            end
        end

    end

end
