import FastGaussQuadrature
import QuadratureRules: shift_nodes, unshift_nodes

@testset "$(rpad("Gauß-Legendre",80))" begin

    for s in 1:10
        @test GaussLegendreQuadrature(s) ≈  GaussLegendreQuadrature(s; fast=true)
        @test GaussLegendreQuadrature(s) == GaussLegendreQuadrature(Float64, s)

        @test sum(weights(GaussLegendreQuadrature(s))) ≈ 1

        @test gauss_legendre_nodes(Float64, s) == nodes(GaussLegendreQuadrature(s))
        @test gauss_legendre_nodes(s) == gauss_legendre_nodes(Float64, s)
        @test gauss_legendre_points(s) == gauss_legendre_points(Float64, s)

        # the points are primary and the nodes are derived from them, so the two agree
        # exactly at equal working precision and up to rounding across precisions
        @test shift_nodes(gauss_legendre_points(BigFloat, s)) == gauss_legendre_nodes(BigFloat, s)
        @test shift_nodes(gauss_legendre_points(Float64, s))  ≈  gauss_legendre_nodes(Float64, s)
        @test unshift_nodes(gauss_legendre_nodes(Float64, s)) ≈ gauss_legendre_points(Float64, s)
        @test gauss_legendre_points(Float64, s) ≈ FastGaussQuadrature.gausslegendre(s)[1]

        @test eltype(gauss_legendre_nodes(BigFloat, s)) == BigFloat
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
