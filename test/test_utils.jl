import QuadratureRules: @big, _big
import QuadratureRules: scale_weights, shift_nodes, shift!, unshift!

@testset "$(rpad("Utility Functions",80))" begin
    
    @test _big(1)    == BigFloat(1)
    @test _big(1//1) == BigFloat(1)
    @test _big(1.0)  == BigFloat(1)
    @test _big("1")  == BigFloat(1)
    @test _big(:(1)) == BigFloat(1)

    x = @big [1  2.0  1//3]
    y = big.([1  2.0  1//3])
    z = [big(1)  big(2.0)  big(1//3)]

    @test x != y
    @test x == z


    b = rand(5); b̃ = copy(b)
    c = rand(5); c̃ = copy(c)

    shift!(b̃, c̃)

    @test b̃ == scale_weights(b)
    @test c̃ == shift_nodes(c)

    unshift!(b̃, c̃)

    @test b̃ == b
    @test c̃ == c

end
