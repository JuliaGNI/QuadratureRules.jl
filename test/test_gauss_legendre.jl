@testset "$(rpad("Gauß-Legendre",80))" begin

    for s in 1:10
        @test GaussLegendreQuadrature(s) ≈  GaussLegendreQuadrature(s; fast=true)
        @test GaussLegendreQuadrature(s) == GaussLegendreQuadrature(Float64, s)

        @test sum(weights(GaussLegendreQuadrature(s))) ≈ 1
    end

    # The default constructor computes nodes and weights in arbitrary precision.
    # A Gauß rule with s nodes integrates polynomials up to degree 2s-1 exactly,
    # so the residual has to be at BigFloat and not at Float64 precision.
    for s in 1:10
        q = GaussLegendreQuadrature(BigFloat, s)
        c = nodes(q)
        b = weights(q)

        @test eltype(c) == BigFloat
        @test eltype(b) == BigFloat

        @test all(0 ≤ cᵢ ≤ 1 for cᵢ in c)
        @test issorted(c, lt = <)

        for k in 0:2s-1
            @test sum(b .* c.^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
        end
    end

end
