import FastGaussQuadrature
import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes

@testset "$(rpad("Radau-Legendre",80))" begin

    # Radau rules, unlike Lobatto rules, are defined for a single node, where they
    # degenerate into the corresponding Riemann sum
    @test RadauLegendreQuadrature(1, :left) == RiemannQuadratureLeft()
    @test RadauLegendreQuadrature(1, :right) == RiemannQuadratureRight()

    # keywords are forwarded, so an unsupported one is an error rather than ignored
    @test_throws MethodError radau_legendre_nodes(4, :left; nosuchkeyword = 1, interval = SymmetricInterval())
    @test_throws MethodError radau_legendre_nodes(4, :left; nosuchkeyword = 1)
    @test_throws MethodError radau_legendre_weights(4, :left; nosuchkeyword = 1)
    @test_throws MethodError radau_legendre_weights(4, :left; nosuchkeyword = 1, interval = SymmetricInterval())
    @test_throws MethodError RadauLegendreQuadrature(4, :left; nosuchkeyword = 1)

    for endpoint in (:left, :right)
        for s in 1:10
            @test RadauLegendreQuadrature(s, endpoint) ≈
                  RadauLegendreQuadrature(s, endpoint; fast = true)
            @test RadauLegendreQuadrature(s, endpoint) ==
                  RadauLegendreQuadrature(Float64, s, Val(endpoint))

            @test sum(weights(RadauLegendreQuadrature(s, endpoint))) ≈ 1

            @test radau_legendre_nodes(Float64, s, Val(endpoint)) ==
                  nodes(RadauLegendreQuadrature(s, endpoint))
            @test radau_legendre_weights(Float64, s, Val(endpoint)) ==
                  weights(RadauLegendreQuadrature(s, endpoint))

            @test radau_legendre_nodes(s, endpoint; interval = SymmetricInterval()) ==
                  radau_legendre_nodes(Float64, s, Val(endpoint); interval = SymmetricInterval())
            @test radau_legendre_nodes(s, endpoint) ==
                  radau_legendre_nodes(Float64, s, Val(endpoint))
            @test radau_legendre_weights(s, endpoint; interval = SymmetricInterval()) ==
                  radau_legendre_weights(Float64, s, Val(endpoint); interval = SymmetricInterval())
            @test radau_legendre_weights(s, endpoint) ==
                  radau_legendre_weights(Float64, s, Val(endpoint))

            # the symmetric nodes are primary and the unit-interval ones are derived from them,
            # so the two agree exactly at equal working precision and up to rounding across precisions
            @test shift_nodes(radau_legendre_nodes(
                BigFloat, s, Val(endpoint); interval = SymmetricInterval())) ==
                  radau_legendre_nodes(BigFloat, s, Val(endpoint))
            @test shift_nodes(radau_legendre_nodes(Float64, s, Val(endpoint); interval = SymmetricInterval())) ≈
                  radau_legendre_nodes(Float64, s, Val(endpoint))
            @test unshift_nodes(radau_legendre_nodes(Float64, s, Val(endpoint))) ≈
                  radau_legendre_nodes(Float64, s, Val(endpoint); interval = SymmetricInterval())

            @test scale_weights(radau_legendre_weights(
                BigFloat, s, Val(endpoint); interval = SymmetricInterval())) ==
                  radau_legendre_weights(BigFloat, s, Val(endpoint))
            @test unscale_weights(radau_legendre_weights(Float64, s, Val(endpoint))) ≈
                  radau_legendre_weights(Float64, s, Val(endpoint); interval = SymmetricInterval())

            @test sum(radau_legendre_weights(Float64, s, Val(endpoint); interval = SymmetricInterval())) ≈
                  2
            @test all(wᵢ > 0
            for wᵢ in radau_legendre_weights(Float64, s, Val(endpoint); interval = SymmetricInterval()))

            @test eltype(radau_legendre_nodes(BigFloat, s, Val(endpoint))) == BigFloat
            @test eltype(radau_legendre_weights(BigFloat, s, Val(endpoint))) == BigFloat
        end
    end

    for s in 1:10
        # the left variant is the classical Gauß-Radau rule, which fixes the first node
        @test radau_legendre_nodes(Float64, s, Val(:left); interval = SymmetricInterval()) ≈
              FastGaussQuadrature.gaussradau(s)[1]
        @test radau_legendre_weights(Float64, s, Val(:left); interval = SymmetricInterval()) ≈
              FastGaussQuadrature.gaussradau(s)[2]

        # the right variant is the reflection of the left one, exactly so
        @test radau_legendre_nodes(BigFloat, s, Val(:right); interval = SymmetricInterval()) ==
              -reverse(radau_legendre_nodes(BigFloat, s, Val(:left); interval = SymmetricInterval()))
        @test radau_legendre_weights(BigFloat, s, Val(:right); interval = SymmetricInterval()) ==
              reverse(radau_legendre_weights(BigFloat, s, Val(:left); interval = SymmetricInterval()))

        # the prescribed endpoint is pinned exactly, and shift_nodes maps it to exactly 0 or 1
        @test radau_legendre_nodes(Float64, s, Val(:left); interval = SymmetricInterval())[begin] ==
              -1
        @test radau_legendre_nodes(Float64, s, Val(:right); interval = SymmetricInterval())[end] ==
              +1
        @test radau_legendre_nodes(Float64, s, Val(:left))[begin] == 0
        @test radau_legendre_nodes(Float64, s, Val(:right))[end] == 1

        # only one endpoint is prescribed, the other one is free and interior
        @test radau_legendre_nodes(Float64, s, Val(:left); interval = SymmetricInterval())[end] <
              +1
        @test radau_legendre_nodes(Float64, s, Val(:right); interval = SymmetricInterval())[begin] >
              -1
    end

    # The default constructor computes nodes and weights in arbitrary precision.
    # A Radau rule with s nodes integrates polynomials up to degree 2s-2 exactly,
    # so the residual has to be at BigFloat and not at Float64 precision.
    for endpoint in (:left, :right)
        for s in 1:10
            q = RadauLegendreQuadrature(BigFloat, s, Val(endpoint))
            c = nodes(q)
            b = weights(q)

            @test eltype(c) == BigFloat
            @test eltype(b) == BigFloat

            @test all(0 ≤ cᵢ ≤ 1 for cᵢ in c)
            @test issorted(c, lt = <)

            for k in 0:(2s - 2)
                @test sum(b .* c .^ k) ≈ 1 / BigFloat(k+1) atol=1E-60
            end
        end
    end
end
