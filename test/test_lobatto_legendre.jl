import FastGaussQuadrature
import QuadratureRules: shift_nodes, unshift_nodes

@testset "$(rpad("Lobatto-Legendre",80))" begin

    @test_throws ErrorException LobattoLegendreQuadrature(1)
    @test_throws ErrorException lobatto_legendre_nodes(1)
    @test_throws ErrorException lobatto_legendre_points(1)

    for s in 2:10
        @test LobattoLegendreQuadrature(s) ≈  LobattoLegendreQuadrature(s; fast=true)
        @test LobattoLegendreQuadrature(s) == LobattoLegendreQuadrature(Float64, s)

        @test sum(weights(LobattoLegendreQuadrature(s))) ≈ 1

        @test lobatto_legendre_nodes(Float64, s) == nodes(LobattoLegendreQuadrature(s))
        @test lobatto_legendre_nodes(s) == lobatto_legendre_nodes(Float64, s)
        @test lobatto_legendre_points(s) == lobatto_legendre_points(Float64, s)

        # the points are primary and the nodes are derived from them, so the two agree
        # exactly at equal working precision and up to rounding across precisions
        @test shift_nodes(lobatto_legendre_points(BigFloat, s)) == lobatto_legendre_nodes(BigFloat, s)
        @test shift_nodes(lobatto_legendre_points(Float64, s))  ≈  lobatto_legendre_nodes(Float64, s)
        @test unshift_nodes(lobatto_legendre_nodes(Float64, s)) ≈ lobatto_legendre_points(Float64, s)
        @test lobatto_legendre_points(Float64, s) ≈ FastGaussQuadrature.gausslobatto(s)[1]

        # the endpoints stay exact under the shift to [-1,+1]
        @test lobatto_legendre_points(Float64, s)[begin] == -1
        @test lobatto_legendre_points(Float64, s)[end]   == +1

        @test lobatto_legendre_nodes(Float64, s; fast=true) ≈ lobatto_legendre_nodes(Float64, s)
        @test eltype(lobatto_legendre_nodes(BigFloat, s)) == BigFloat
    end

    # The default constructor computes nodes and weights in arbitrary precision.
    # A Lobatto rule with s nodes integrates polynomials up to degree 2s-3
    # exactly, so the residual has to be at BigFloat and not at Float64
    # precision. The first and last node are the endpoints of the interval.
    for s in 2:10
        q = LobattoLegendreQuadrature(BigFloat, s)
        c = nodes(q)
        b = weights(q)

        @test eltype(c) == BigFloat
        @test eltype(b) == BigFloat

        @test c[begin] == 0
        @test c[end]   == 1
        @test issorted(c, lt = <)

        for k in 0:2s-3
            @test sum(b .* c.^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
        end
    end

end
