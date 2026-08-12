import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes, shift!, unshift!
import QuadratureRules: _default_arithmetic

# A stand-in for an element type from outside the numeric tower, such as that of a computer
# algebra system, of which the default arithmetic knows nothing but that it is not numeric.
struct NotANumber end

@testset "$(rpad("Utility Functions",80))" begin

    # The numeric tower is computed in BigFloat and rounded to the element type at the very
    # end. Rational and Complex are counted as numeric, not as exact: a quadrature node is an
    # algebraic number that no rational represents, and the exact root find needs `eigvals`.
    @test _default_arithmetic(Float16)          == BigFloat
    @test _default_arithmetic(Float32)          == BigFloat
    @test _default_arithmetic(Float64)          == BigFloat
    @test _default_arithmetic(BigFloat)         == BigFloat
    @test _default_arithmetic(Int)              == BigFloat
    @test _default_arithmetic(BigInt)           == BigFloat
    @test _default_arithmetic(Rational{Int})    == BigFloat
    @test _default_arithmetic(Rational{BigInt}) == BigFloat
    @test _default_arithmetic(ComplexF64)       == BigFloat

    # everything else is taken to do its own arithmetic exactly and is computed in itself
    @test _default_arithmetic(NotANumber) == NotANumber

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
