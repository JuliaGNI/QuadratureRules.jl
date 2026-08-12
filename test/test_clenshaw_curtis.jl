import FastTransforms
import QuadratureRules: scale_weights, unscale_weights, shift_nodes, unshift_nodes

@testset "$(rpad("Clenshaw-Curtis",80))" begin

    @test_throws ErrorException ClenshawCurtisQuadrature(1)
    @test_throws ErrorException clenshaw_curtis_weights(1)
    @test_throws ErrorException clenshaw_curtis_weights(1; interval = SymmetricInterval())

    for s in 2:10
        @test ClenshawCurtisQuadrature(s) == ClenshawCurtisQuadrature(Float64, s)

        @test sum(weights(ClenshawCurtisQuadrature(s))) ≈ 1
        @test order(ClenshawCurtisQuadrature(s)) == (isodd(s) ? s+1 : s)

        # an s-node interpolatory rule integrates polynomials up to degree s-1 exactly,
        # and one degree further for odd s; see test_order.jl
        let q = ClenshawCurtisQuadrature(BigFloat, s)
            for k in 0:order(q)-1
                @test sum(weights(q) .* nodes(q).^k) ≈ 1 / BigFloat(k+1)  atol=1E-60
            end
        end


        μ = FastTransforms.chebyshevmoments1(Float64, s)
        b = scale_weights(reverse(FastTransforms.clenshawcurtisweights(μ)))
        c = shift_nodes(reverse(FastTransforms.clenshawcurtisnodes(Float64, s)))

        @test nodes(ClenshawCurtisQuadrature(s)) ≈ c
        @test weights(ClenshawCurtisQuadrature(s)) ≈ b

        @test clenshaw_curtis_nodes(Float64, s) ≈ c

        # both compute the nodes in the working precision IT and round to T, so
        # they agree exactly and not merely approximately
        @test clenshaw_curtis_nodes(Float64, s)  ==  nodes(ClenshawCurtisQuadrature(s))
        @test clenshaw_curtis_nodes(Float32, s)  ==  nodes(ClenshawCurtisQuadrature(Float32, s))
        @test clenshaw_curtis_nodes(Float64, s; IT=Float64) == nodes(ClenshawCurtisQuadrature(Float64, s; IT=Float64))
        @test clenshaw_curtis_nodes(BigFloat, s) == nodes(ClenshawCurtisQuadrature(BigFloat, s))
        @test clenshaw_curtis_nodes(Float64, s; interval = SymmetricInterval()) ≈ reverse(FastTransforms.clenshawcurtisnodes(Float64, s))
        @test unshift_nodes(clenshaw_curtis_nodes(Float64, s)) ≈ clenshaw_curtis_nodes(Float64, s; interval = SymmetricInterval())

        @test clenshaw_curtis_nodes(s; interval = SymmetricInterval()) == clenshaw_curtis_nodes(Float64, s; interval = SymmetricInterval())
        @test clenshaw_curtis_nodes(s)  == clenshaw_curtis_nodes(Float64, s)

        # the same holds for the weights, which are primary on [0,1] here
        @test clenshaw_curtis_weights(Float64, s) ≈ b
        @test clenshaw_curtis_weights(Float64, s)  ==  weights(ClenshawCurtisQuadrature(s))
        @test clenshaw_curtis_weights(Float32, s)  ==  weights(ClenshawCurtisQuadrature(Float32, s))
        @test clenshaw_curtis_weights(Float64, s; IT=Float64) == weights(ClenshawCurtisQuadrature(Float64, s; IT=Float64))
        @test clenshaw_curtis_weights(BigFloat, s) == weights(ClenshawCurtisQuadrature(BigFloat, s))

        @test clenshaw_curtis_weights(s)       == clenshaw_curtis_weights(Float64, s)
        @test clenshaw_curtis_weights(s; interval = SymmetricInterval()) == clenshaw_curtis_weights(Float64, s; interval = SymmetricInterval())

        @test unscale_weights(clenshaw_curtis_weights(BigFloat, s)) == clenshaw_curtis_weights(BigFloat, s; interval = SymmetricInterval())
        @test scale_weights(clenshaw_curtis_weights(Float64, s; interval = SymmetricInterval())) ≈ clenshaw_curtis_weights(Float64, s)
        @test sum(clenshaw_curtis_weights(Float64, s; interval = SymmetricInterval())) ≈ 2

        # all Clenshaw-Curtis weights are positive (Imhof, 1963); this is what
        # guarantees convergence for every continuous integrand
        @test all(w -> w > 0, weights(ClenshawCurtisQuadrature(s)))

        # the IT keyword selects the working precision; all choices must agree
        @test ClenshawCurtisQuadrature(Float64, s; IT=Float64) ≈ ClenshawCurtisQuadrature(s)
        @test LobattoChebyshevQuadrature(Float64, s; IT=Float64) ≈ LobattoChebyshevQuadrature(s)
        @test eltype(ClenshawCurtisQuadrature(Float32, s; IT=Float32)) == Float32
    end

    # Reid's explicit closed form for the weights on [-1,+1] for even N (N = s-1),
    # 18.330 Lecture Notes: Clenshaw-Curtis Quadrature, Section 4. Note that the
    # last term of the sum carries a factor 1/2, which Reid's sample tables omit;
    # without it the rule loses its degree-N exactness for even N.
    function reid_weights(N)
        w = zeros(N + 1)
        w[1] = w[N+1] = 1 / (N^2 - 1)
        for m in 1:N-1
            t = sum(n -> 2 / (1 - 4n^2) * cos(2 * m * n * π / N), 1:(N÷2-1); init = 0.0)
            w[m+1] = 2 / N * (1 + t + cos(m * π) / (1 - N^2))
        end
        return w
    end

    for N in 4:2:16
        # the package normalises to [0,1], Reid to [-1,+1], hence the factor 2;
        # Reid indexes his nodes from t = +1 downwards, hence the reverse
        @test 2 .* reverse(weights(ClenshawCurtisQuadrature(N+1))) ≈ reid_weights(N)
    end

    # Aliasing on the Chebyshev grid: Tₙ₊ₚ and Tₙ₋ₚ are indistinguishable at the
    # nodes, so the rule returns the exact integral of Tₙ₋ₚ when given Tₙ₊ₚ
    # (Trefethen, 2008, Theorem 5.2 and equation 5.3).
    chebT(j, x) = cos(j * acos(clamp(x, -1, 1)))
    unit_to_symmetric(quad, f) = 2 * quad(ξ -> f(2ξ - 1))

    for n in 4:12, p in 0:min(3, n)
        quad = ClenshawCurtisQuadrature(n+1)
        m = n - p
        expected = isodd(m) ? 0.0 : 2 / (1 - m^2)
        @test unit_to_symmetric(quad, x -> chebT(n+p, x)) ≈ expected  atol=1E-12
    end

end
