import GeometricBase

@testset "$(rpad("QuadratureRule",80))" begin
    o = 2
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    quad = QuadratureRule(o, c, b)

    @test order(quad) == o
    @test nodes(quad) == c
    @test weights(quad) == b
    @test nnodes(quad) == length(b)
    @test eachindex(quad) == 1:2
    @test eltype(quad) == Float64

    @test quad(x -> 1) == 1
    @test quad(x -> x) == 1//2
    @test quad(x -> x^2) == 1//2

    quad2 = QuadratureRule{Float64, 2}(o, c, b)
    quadϵ = QuadratureRule{Float64, 2}(o, c .+ 8eps() * rand(2), b .+ 8eps() * rand(2))

    @test hash(quad) == hash(quad2)
    @test hash(quad) != hash(quadϵ)

    @test isequal(quad, quad2)
    @test !isequal(quad, quadϵ)

    @test isapprox(quad, quad2)
    @test isapprox(quad, quadϵ)

    @test quad ≈ quad2
    @test quad == quad2
    @test quad === quad2

    @test quad ≈ quadϵ
    @test quad != quadϵ
    @test quad !== quadϵ
end

@testset "$(rpad("Accessor bindings",80))" begin

    # the accessors must be methods on the GeometricBase generics, not functions of their
    # own, or loading QuadratureRules together with another package of the ecosystem makes
    # the exported names resolve to nothing

    @test QuadratureRules.nnodes === GeometricBase.nnodes
    @test QuadratureRules.nodes === GeometricBase.nodes
    @test QuadratureRules.order === GeometricBase.order
    @test QuadratureRules.weights === GeometricBase.weights
end
