import FastTransforms
import QuadratureRules: shift_nodes

@testset "$(rpad("Lobatto-Chebyshev",80))" begin

    for s in 2:10
        c = shift_nodes(reverse(FastTransforms.chebyshevpoints(Float64, s, Val(2)), 1))

        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(s, Val(2))
        @test LobattoChebyshevQuadrature(s) == ChebyshevQuadrature(Float64, s, Val(2))
        
        @test nodes(LobattoChebyshevQuadrature(s)) ≈ c

        # TODO test weights

    end

end
