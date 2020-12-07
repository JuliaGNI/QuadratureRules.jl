@testset "$(rpad("Gauß-Legendre",80))" begin

    for s in 1:10
        @test GaussLegendreQuadrature(s) ≈  GaussLegendreQuadrature(s; fast=true)
        @test GaussLegendreQuadrature(s) == GaussLegendreQuadrature(Float64, s)

        @test sum(weights(GaussLegendreQuadrature(s))) ≈ 1
    end

end
