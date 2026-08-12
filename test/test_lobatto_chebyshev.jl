import FastTransforms
import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes

@testset "$(rpad("Lobatto-Chebyshev",80))" begin

    @test_throws ErrorException LobattoChebyshevQuadrature(1)
    @test_throws ErrorException chebyshev_nodes(Float64, 1, Val(2); interval = SymmetricInterval())
    @test_throws ErrorException chebyshev_nodes(Float64, 1, Val(2))
    @test_throws ErrorException chebyshev_weights(Float64, 1, Val(2))
    @test_throws ErrorException chebyshev_weights(Float64, 1, Val(2); interval = SymmetricInterval())
    @test_throws ErrorException lobatto_chebyshev_nodes(1; interval = SymmetricInterval())
    @test_throws ErrorException lobatto_chebyshev_nodes(1)
    @test_throws ErrorException lobatto_chebyshev_weights(1)
    @test_throws ErrorException lobatto_chebyshev_weights(1; interval = SymmetricInterval())

    for s in 2:10
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(s, 2)
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(2))
        
        @test sum(weights(LobattoChebyshevQuadrature(s))) ≈ 1

        x = reverse(FastTransforms.chebyshevpoints(Float64, s, Val(2)))
        c = shift_nodes(x)

        @test chebyshev_nodes(Float64, s, Val(2); interval = SymmetricInterval()) ≈ x
        @test chebyshev_nodes(Float64, s, Val(2))  ≈ c
        @test nodes(LobattoChebyshevQuadrature(s)) ≈ c

        @test lobatto_chebyshev_nodes(Float64, s; interval = SymmetricInterval()) == chebyshev_nodes(Float64, s, Val(2); interval = SymmetricInterval())
        @test lobatto_chebyshev_nodes(Float64, s)  == chebyshev_nodes(Float64, s, Val(2))
        @test lobatto_chebyshev_nodes(s; interval = SymmetricInterval()) == lobatto_chebyshev_nodes(Float64, s; interval = SymmetricInterval())
        @test lobatto_chebyshev_nodes(s)  == lobatto_chebyshev_nodes(Float64, s)
        @test chebyshev_nodes(s, 2; interval = SymmetricInterval()) == chebyshev_nodes(Float64, s, Val(2); interval = SymmetricInterval())
        @test chebyshev_nodes(s, 2)  == chebyshev_nodes(Float64, s, Val(2))

        @test unshift_nodes(chebyshev_nodes(Float64, s, Val(2))) ≈ chebyshev_nodes(Float64, s, Val(2); interval = SymmetricInterval())

        @test chebyshev_weights(Float64, s, Val(2)) == weights(LobattoChebyshevQuadrature(s))
        @test lobatto_chebyshev_weights(Float64, s)       == chebyshev_weights(Float64, s, Val(2))
        @test lobatto_chebyshev_weights(Float64, s; interval = SymmetricInterval()) == chebyshev_weights(Float64, s, Val(2); interval = SymmetricInterval())
        @test lobatto_chebyshev_weights(s)       == lobatto_chebyshev_weights(Float64, s)
        @test lobatto_chebyshev_weights(s; interval = SymmetricInterval()) == lobatto_chebyshev_weights(Float64, s; interval = SymmetricInterval())
        @test chebyshev_weights(s, 2)       == chebyshev_weights(Float64, s, Val(2))
        @test chebyshev_weights(s, 2; interval = SymmetricInterval()) == chebyshev_weights(Float64, s, Val(2); interval = SymmetricInterval())

        # here the weights for [0,1] are primary and those for [-1,+1] are derived
        @test unscale_weights(chebyshev_weights(BigFloat, s, Val(2))) == chebyshev_weights(BigFloat, s, Val(2); interval = SymmetricInterval())
        @test scale_weights(chebyshev_weights(Float64, s, Val(2); interval = SymmetricInterval())) ≈ chebyshev_weights(Float64, s, Val(2))
        @test sum(chebyshev_weights(Float64, s, Val(2); interval = SymmetricInterval())) ≈ 2

        # the Chebyshev points of the second kind are exactly the Clenshaw-Curtis
        # nodes, so the interpolatory rule on them is the Clenshaw-Curtis rule
        @test LobattoChebyshevQuadrature(s) == ClenshawCurtisQuadrature(s)
        @test lobatto_chebyshev_weights(BigFloat, s) == clenshaw_curtis_weights(BigFloat, s)
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
