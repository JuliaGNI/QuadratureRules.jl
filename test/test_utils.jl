import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes, shift!, unshift!

@testset "$(rpad("Utility Functions",80))" begin

    b = rand(5); b̃ = copy(b)
    c = rand(5); c̃ = copy(c)

    shift!(b̃, c̃)

    @test b̃ ≈ scale_weights(b)  atol=1E-14
    @test c̃ ≈ shift_nodes(c)    atol=1E-14

    unshift!(b̃, c̃)

    @test b̃ ≈ b  atol=1E-14
    @test c̃ ≈ c  atol=1E-14

    @test unshift_nodes(shift_nodes(c)) ≈ c  atol=1E-14
    @test shift_nodes(unshift_nodes(c)) ≈ c  atol=1E-14
    @test unshift_nodes([0.0, 1.0]) == [-1.0, 1.0]

    @test unscale_weights(scale_weights(b)) ≈ b  atol=1E-14
    @test scale_weights(unscale_weights(b)) ≈ b  atol=1E-14
    @test unscale_weights([0.5, 0.5]) == [1.0, 1.0]

end
