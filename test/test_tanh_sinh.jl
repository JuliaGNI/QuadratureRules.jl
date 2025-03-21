@testset "$(rpad("Tanh-Sinh",80))" begin

    @test_throws ErrorException TanhSinhQuadrature(0)

    q = TanhSinhQuadrature(1)
    @test sum(weights(q)) ≈ 1  atol=1e-5

    q = TanhSinhQuadrature(2)
    @test sum(weights(q)) ≈ 1  atol=1e-13

    q = TanhSinhQuadrature(3)
    @test sum(weights(q)) ≈ 1  atol=1e-15

end
