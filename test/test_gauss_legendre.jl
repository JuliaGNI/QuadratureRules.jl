import FastGaussQuadrature
import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes

@testset "$(rpad("Gauß-Legendre",80))" begin

    for s in 1:10
        @test GaussLegendreQuadrature(s) ≈  GaussLegendreQuadrature(s; fast=true)
        @test GaussLegendreQuadrature(s) == GaussLegendreQuadrature(Float64, s)

        @test sum(weights(GaussLegendreQuadrature(s))) ≈ 1

        @test gauss_legendre_nodes(Float64, s) == nodes(GaussLegendreQuadrature(s))
        @test gauss_legendre_nodes(s) == gauss_legendre_nodes(Float64, s)
        @test symmetric_gauss_legendre_nodes(s) == symmetric_gauss_legendre_nodes(Float64, s)

        @test gauss_legendre_weights(Float64, s) == weights(GaussLegendreQuadrature(s))
        @test gauss_legendre_weights(s) == gauss_legendre_weights(Float64, s)
        @test symmetric_gauss_legendre_weights(s) == symmetric_gauss_legendre_weights(Float64, s)

        # the symmetric nodes are primary and the unit-interval ones are derived from them,
        # so the two agree exactly at equal working precision and up to rounding across precisions
        @test shift_nodes(symmetric_gauss_legendre_nodes(BigFloat, s)) == gauss_legendre_nodes(BigFloat, s)
        @test shift_nodes(symmetric_gauss_legendre_nodes(Float64, s))  ≈  gauss_legendre_nodes(Float64, s)
        @test unshift_nodes(gauss_legendre_nodes(Float64, s)) ≈ symmetric_gauss_legendre_nodes(Float64, s)
        @test symmetric_gauss_legendre_nodes(Float64, s) ≈ FastGaussQuadrature.gausslegendre(s)[1]

        # likewise the weights for [-1,+1] are primary and those for [0,1] are derived
        @test scale_weights(symmetric_gauss_legendre_weights(BigFloat, s)) == gauss_legendre_weights(BigFloat, s)
        @test unscale_weights(gauss_legendre_weights(Float64, s))      ≈  symmetric_gauss_legendre_weights(Float64, s)
        @test symmetric_gauss_legendre_weights(Float64, s) ≈ FastGaussQuadrature.gausslegendre(s)[2]

        @test sum(symmetric_gauss_legendre_weights(Float64, s)) ≈ 2
        @test all(wᵢ > 0 for wᵢ in symmetric_gauss_legendre_weights(Float64, s))

        @test eltype(gauss_legendre_nodes(BigFloat, s)) == BigFloat
        @test eltype(gauss_legendre_weights(BigFloat, s)) == BigFloat
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

    # Element types from the numeric tower other than the floating point ones are computed in
    # BigFloat like those, and not in themselves as a symbolic type would be: their arithmetic
    # is of no help for an algebraic number, and the exact root find needs `eigvals`, which
    # they do not have.
    for s in 1:10
        @test symmetric_gauss_legendre_nodes(ComplexF64, s) == ComplexF64.(symmetric_gauss_legendre_nodes(BigFloat, s))
        @test gauss_legendre_weights(ComplexF64, s) == ComplexF64.(gauss_legendre_weights(BigFloat, s))

        @test symmetric_gauss_legendre_nodes(Rational{BigInt}, s) == Rational{BigInt}.(symmetric_gauss_legendre_nodes(BigFloat, s))
        @test gauss_legendre_weights(Rational{BigInt}, s) == Rational{BigInt}.(gauss_legendre_weights(BigFloat, s))
    end

end
