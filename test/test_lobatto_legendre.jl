@testset "$(rpad("Lobatto-Legendre",80))" begin

    @test_throws ErrorException LobattoLegendreQuadrature(1)

    for s in 2:10
        @test LobattoLegendreQuadrature(s) ≈  LobattoLegendreQuadrature(s; fast=true)
        @test LobattoLegendreQuadrature(s) == LobattoLegendreQuadrature(Float64, s)

        @test sum(weights(LobattoLegendreQuadrature(s))) ≈ 1
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
