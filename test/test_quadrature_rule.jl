@testset "$(rpad("QuadratureRule",80))" begin

    o = 2
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    quad = QuadratureRule(o,c,b)

    @test order(quad) == o
    @test nodes(quad) == c
    @test weights(quad) == b
    @test nnodes(quad) == length(b)
    @test eachindex(quad) == 1:2

    @test quad(x -> 1)   == 1
    @test quad(x -> x)   == 1//2
    @test quad(x -> x^2) == 1//2


    quad2 = QuadratureRule{Float64,2}(o, c, b)
    quadϵ = QuadratureRule{Float64,2}(o, c .+ 8eps() * rand(2), b .+ 8eps() * rand(2))

    @test  isequal(quad, quad2)
    @test !isequal(quad, quadϵ)

    @test isapprox(quad, quad2)
    @test isapprox(quad, quadϵ)

    @test quad ≈   quad2
    @test quad ==  quad2
    @test quad === quad2

    @test quad ≈   quadϵ
    @test quad !=  quadϵ
    @test quad !== quadϵ

end
