@testset "$(rpad("Lobatto-Legendre",80))" begin

    @test_throws ErrorException LobattoLegendreQuadrature(1)

    for s in 2:10
        @test LobattoLegendreQuadrature(s) ≈  LobattoLegendreQuadrature(s; fast=true)
        @test LobattoLegendreQuadrature(s) == LobattoLegendreQuadrature(Float64, s)

        @test sum(weights(LobattoLegendreQuadrature(s))) ≈ 1
    end

end
