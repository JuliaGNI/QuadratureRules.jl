
@testset "$(rpad("Order",80))" begin

    # The order p of a rule is sharp: it integrates every polynomial of degree ≤ p-1
    # exactly and fails for some polynomial of degree p. Since integration is linear it
    # suffices to test the monomials, and since the nodes lie in [0,1] the k-th moment
    # is 1/(k+1). Everything is evaluated in BigFloat so that the two assertions below
    # separate the property of the rule from the round-off of the arithmetic.
    moment(quad, k) = sum(weights(quad) .* nodes(quad).^k)

    exact_to(quad, k)  = isapprox(moment(quad, k), 1 / BigFloat(k+1), atol = 1E-60)
    exact_upto(quad)   = all(exact_to(quad, k) for k in 0:order(quad)-1)
    sharp(quad)        = !exact_to(quad, order(quad))

    # the second assertion is the one that pins the order down. Without it an
    # understated order passes unnoticed, which is how the Chebyshev rules came to
    # report s where they are in fact exact to degree s for odd s.
    function test_order(quad)
        @test exact_upto(quad)
        @test sharp(quad)
    end


    @testset "$(rpad("tabulated rules",60))" begin
        for quadrature in (RiemannQuadratureLeft, RiemannQuadratureRight,
                           MidpointQuadrature, TrapezoidalQuadrature)
            test_order(quadrature(BigFloat))
        end
    end


    @testset "$(rpad("generated rules",60))" begin
        for s in 2:12
            test_order(GaussLegendreQuadrature(BigFloat, s))
            test_order(LobattoLegendreQuadrature(BigFloat, s))
            test_order(ClenshawCurtisQuadrature(BigFloat, s))
            test_order(GaussChebyshevQuadrature(BigFloat, s))
            test_order(LobattoChebyshevQuadrature(BigFloat, s))
        end

        # Gauss-Chebyshev is the only generated rule defined for a single node
        test_order(GaussChebyshevQuadrature(BigFloat, 1))
    end


    @testset "$(rpad("reported orders",60))" begin
        for s in 2:12
            @test order(GaussLegendreQuadrature(s))     == 2s
            @test order(LobattoLegendreQuadrature(s))   == 2s-2

            # the Chebyshev rules pick up one degree for an odd number of nodes,
            # because the monomial of degree s they would otherwise fail on is odd
            # about the midpoint of the interval
            @test order(ClenshawCurtisQuadrature(s))    == (isodd(s) ? s+1 : s)
            @test order(GaussChebyshevQuadrature(s))    == (isodd(s) ? s+1 : s)
            @test order(LobattoChebyshevQuadrature(s))  == (isodd(s) ? s+1 : s)
        end

        @test order(GaussChebyshevQuadrature(1)) == 2
    end


    # Tanh-sinh is the one rule in the package without a degree of exactness. It is the
    # trapezoidal rule after a change of variables, not an interpolatory rule, and its
    # truncated sum reproduces no polynomial exactly — not even the constant, so `moment`
    # already fails at k = 0. Its order is therefore 0, and it is absent from the two loops
    # above: `test_order` would assert exactness up to degree -1, and `exact_upto` over the
    # empty range `0:-1` is vacuously true, which would assert nothing at all.
    @testset "$(rpad("tanh-sinh has no order",60))" begin
        for n in 1:5
            quad = TanhSinhQuadrature(BigFloat, n)

            @test order(quad) == 0
            @test moment(quad, 0) != 1        # never exact, however close it gets
            @test moment(quad, 0) ≈ 1  atol = 1E-5

            # From level 4 the truncation error falls below the 1E-60 tolerance that
            # `exact_to` uses throughout this file, so the moment machinery can no longer
            # tell the rule from an exact one. That is another reason why the accuracy of
            # tanh-sinh is tested by its convergence rate in test_tanh_sinh.jl instead.
            n ≤ 3 && @test !exact_to(quad, 0)
        end
    end


    # Because the order is sharp, a rule is determined by its nodes and weights alone,
    # so rules that coincide compare equal even across families. An interpolatory rule
    # is uniquely determined by its nodes, and these node sets agree exactly.
    @testset "$(rpad("rules that coincide",60))" begin
        @test ClenshawCurtisQuadrature(2)   == TrapezoidalQuadrature()
        @test ClenshawCurtisQuadrature(2)   == LobattoLegendreQuadrature(2)
        @test ClenshawCurtisQuadrature(3)   == LobattoLegendreQuadrature(3)   # Simpson's
        @test LobattoChebyshevQuadrature(3) == LobattoLegendreQuadrature(3)

        @test GaussChebyshevQuadrature(1)   == MidpointQuadrature()
        @test GaussChebyshevQuadrature(1)   == GaussLegendreQuadrature(1)
        @test MidpointQuadrature()          == GaussLegendreQuadrature(1)

        for T in (Float32, Float64)
            @test isequal(ClenshawCurtisQuadrature(T, 3), LobattoLegendreQuadrature(T, 3))
            @test isequal(GaussChebyshevQuadrature(T, 1), MidpointQuadrature(T))
        end

        # At T = BigFloat the working precision IT equals the target precision, so the
        # weights are no longer correctly rounded and the two constructions of Simpson's
        # rule need not agree bit for bit: the Lobatto closed form 2/(s(s-1)Pₛ₋₁²) lands
        # on 1/6 exactly, whereas the Clenshaw-Curtis cosine sum is one ulp below it.
        # This is the IT = T degradation that test_precision.jl quantifies, and it is why
        # BigFloat is the default working precision for the narrower element types.
        let cc = ClenshawCurtisQuadrature(BigFloat, 3), ll = LobattoLegendreQuadrature(BigFloat, 3)
            @test order(cc) == order(ll)
            @test nodes(cc) == nodes(ll)
            @test weights(cc) ≈ weights(ll)  atol = 2eps(big(1)/6)
            @test weights(ll) == [big(1)/6, big(2)/3, big(1)/6]
        end
    end

end
