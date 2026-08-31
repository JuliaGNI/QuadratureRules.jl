import QuadratureRules: scale_weights, unscale_weights, unshift_nodes

@testset "$(rpad("Tanh-Sinh",80))" begin
    @test_throws ErrorException TanhSinhQuadrature(0)
    @test_throws ErrorException TanhSinhQuadrature(-1)
    @test_throws ErrorException tanh_sinh_weights(0)
    @test_throws ErrorException tanh_sinh_weights(0; interval = SymmetricInterval())

    for n in 1:6
        @test TanhSinhQuadrature(n) == TanhSinhQuadrature(Float64, n)
        @test tanh_sinh_nodes(n) == tanh_sinh_nodes(Float64, n)
        @test tanh_sinh_nodes(n; interval = SymmetricInterval()) ==
              tanh_sinh_nodes(Float64, n; interval = SymmetricInterval())
        @test tanh_sinh_weights(n) == tanh_sinh_weights(Float64, n)
        @test tanh_sinh_weights(n; interval = SymmetricInterval()) ==
              tanh_sinh_weights(Float64, n; interval = SymmetricInterval())
    end

    # The level n is a step size h = 2⁻ⁿ, not a node count; the latter follows from where
    # the transformed integrand drops below what T can resolve and grows like 2ⁿ.
    @test [nnodes(TanhSinhQuadrature(n)) for n in 1:5] == [13, 25, 51, 101, 203]

    @testset "$(rpad("structure of the rule",60))" begin
        for T in (Float32, Float64, BigFloat), n in 1:6

            quad = TanhSinhQuadrature(T, n)
            c, b = nodes(quad), weights(quad)

            @test eltype(quad) == T
            @test nnodes(quad) == length(c) == length(b)
            @test isodd(nnodes(quad))                    # a centre node plus symmetric pairs

            # `lt = <` would not catch duplicates, which is precisely the failure mode of a
            # rule whose outer nodes crowd together at the resolution limit of T
            @test issorted(c, lt = <=)
            @test allunique(c)

            # tanh-sinh maps the endpoints to infinity, so no node may land on one: this is
            # what allows an integrand singular at 0 or 1 to be passed in directly
            @test all(0 < cᵢ < 1 for cᵢ in c)
            @test !any(iszero, c)
            @test !any(isone, c)
            @test isfinite(quad(x -> 1 / sqrt(x * (1 - x))))

            @test all(bᵢ > 0 for bᵢ in b)

            # The weights sum to one only up to the truncation error, which is about 10⁻⁶ at
            # level 1 and squares with every further level until eps(T) is reached. See the
            # convergence testset below, which pins that rate down.
            @test sum(b) ≈ one(T) atol = max(4 * eps(T), T(10)^(-5 * 2^(n-1)))

            # nodes and weights are assembled in pairs about the centre node 1/2, so the
            # weights are symmetric bit for bit and the centre node is exactly 1/2
            @test b == reverse(b)
            @test c ≈ reverse(one(T) .- c)
            @test c[cld(nnodes(quad), 2)] == T(1) / 2

            # the nodes on [-1,+1] carry one bit less next to the endpoints than those
            # on [0,1]; the truncation accounts for that, so they too remain distinct
            x = tanh_sinh_nodes(T, n; interval = SymmetricInterval())
            @test issorted(x, lt = <=)
            @test allunique(x)
            @test all(-1 < xᵢ < 1 for xᵢ in x)
            @test x ≈ -reverse(x)

            # both accessors truncate exactly as the rule does and round only once
            @test tanh_sinh_nodes(T, n) == c
            @test unshift_nodes(tanh_sinh_nodes(T, n)) ≈ x
            @test length(x) == length(c)

            # the same holds for the weights, which are primary on [0,1] here, so that the
            # weights for [-1,+1] must be derived from them and truncate alike
            w = tanh_sinh_weights(T, n; interval = SymmetricInterval())
            @test tanh_sinh_weights(T, n) == b
            @test length(w) == length(b)
            @test w == reverse(w)
            @test all(wᵢ > 0 for wᵢ in w)
            @test unscale_weights(tanh_sinh_weights(T, n)) == w
            @test scale_weights(w) ≈ b
            @test sum(w) ≈ 2one(T) atol = 2 * max(4 * eps(T), T(10)^(-5 * 2^(n-1)))
        end
    end

    # Level n uses h = 2⁻ⁿ, so its grid in t is contained in that of level n+1 and hence,
    # the truncation criterion being a threshold in t shared by all levels, so are its nodes.
    # This holds while the rule is still resolving; beyond that the pairwise merge may fold
    # the outermost nodes of two levels together, which is where Float32 gives out at n = 5.
    @testset "$(rpad("nested levels",60))" begin
        for (T, N) in ((Float32, 4), (Float64, 5), (BigFloat, 6)), n in 1:N

            @test issubset(nodes(TanhSinhQuadrature(T, n)), nodes(TanhSinhQuadrature(T, n+1)))
        end
    end

    # Unlike every other rule in this package, tanh-sinh has no degree of exactness: it does
    # not integrate even a constant exactly, the weight sum differing from one by the
    # truncation error. Its reported order is therefore 0, and it is absent from the RULES
    # lists of test_order.jl and test_precision.jl, which assume a positive polynomial order.
    @testset "$(rpad("no polynomial exactness",60))" begin
        for n in 1:6
            @test order(TanhSinhQuadrature(n)) == 0
            @test order(TanhSinhQuadrature(BigFloat, n)) == 0

            # the k = 0 moment, i.e. the integral of the constant 1, is not reproduced.
            # Only up to level 4: beyond it the truncation error falls to a few ulps of 1,
            # so that the inequality would rest on the round-off of the sum rather than on
            # the rule — as it does already in Float32 at level 4 and Float64 at level 6,
            # where the weights do sum to exactly one.
            n ≤ 4 && @test sum(weights(TanhSinhQuadrature(BigFloat, n))) != 1
        end
    end

    # Halving the step size roughly doubles the number of correct digits, until the precision
    # of T is reached. The errors below are pure discretisation errors, independent of the
    # precision: 3.4e-6, 3.7e-14, 1.0e-30, 2.1e-64 for the weight sum at levels 1 to 4.
    @testset "$(rpad("convergence",60))" begin
        INTEGRANDS = ((x -> one(x), big(1)),
            (log, big(-1)),
            (exp, exp(big(1)) - 1),
            (x -> sqrt(x) * log(1/x), big(4)//9))

        for (f, exact) in INTEGRANDS
            err = [abs(TanhSinhQuadrature(BigFloat, n)(f) - exact) for n in 1:4]

            for n in 1:3
                @test err[n + 1] < 1000 * err[n]^2
            end

            @test err[4] < 1E-50
            # and the working precision is reached one level later
            @test abs(TanhSinhQuadrature(BigFloat, 5)(f) - exact) < 1E-74
        end
    end

    # An integrand with a mild endpoint singularity is integrated to full precision.
    @testset "$(rpad("endpoint singularities",60))" begin
        for n in 3:6
            quad = TanhSinhQuadrature(n)
            @test quad(log) ≈ -1 atol=1E-14
            @test quad(x -> sqrt(x) * log(1/x)) ≈ 4/9 atol=1E-14
            @test quad(x -> log(x) * log(1-x)) ≈ 2 - π^2/6 atol=1E-14
        end

        # An algebraic singularity is a different matter. The outermost node cannot approach
        # the endpoint closer than about eps(T), so for an integrand behaving like x^(-1/2)
        # the neglected tail is of size sqrt(eps(T)) — and no further level improves on it.
        # Both bounds are asserted, so that the saturation is pinned down rather than merely
        # tolerated: this is the property that motivates computing in arbitrary precision.
        for T in (Float32, Float64, BigFloat), n in 4:6

            quad = TanhSinhQuadrature(T, n)

            for (f, exact) in ((x -> 1/sqrt(x), T == BigFloat ? big(2) : 2.0),
                (x -> 1/sqrt(x * (1 - x)), T == BigFloat ? big(π) : Float64(π)))
                err = abs(quad(f) - exact)
                @test err ≤ 10 * sqrt(eps(T))
                @test err ≥ 0.1 * sqrt(eps(T))
            end
        end

        # ... and it does improve with the precision of T, which is the point
        @test abs(TanhSinhQuadrature(Float32, 5)(x -> 1/sqrt(x)) - 2) > 1E-5
        @test abs(TanhSinhQuadrature(Float64, 5)(x -> 1/sqrt(x)) - 2) < 1E-7
        @test abs(TanhSinhQuadrature(BigFloat, 5)(x -> 1/sqrt(x)) - big(2)) < 1E-38

        setprecision(BigFloat, 512) do
            quad = TanhSinhQuadrature(BigFloat, 6)
            @test abs(quad(x -> 1/sqrt(x)) - big(2)) < 1E-77
            @test abs(quad(log) + 1) < 1E-150
        end
    end

    # the IT keyword selects the working precision; all choices must agree
    @testset "$(rpad("working precision",60))" begin
        for n in 1:5
            @test TanhSinhQuadrature(Float64, n; IT = Float64) ≈ TanhSinhQuadrature(n)
            @test tanh_sinh_nodes(Float64, n; IT = Float64) ≈ tanh_sinh_nodes(n)
            @test tanh_sinh_nodes(Float64, n; IT = Float64, interval = SymmetricInterval()) ≈
                  tanh_sinh_nodes(n; interval = SymmetricInterval())
            @test tanh_sinh_weights(Float64, n; IT = Float64) ≈ tanh_sinh_weights(n)
            @test tanh_sinh_weights(Float64, n; IT = Float64, interval = SymmetricInterval()) ≈
                  tanh_sinh_weights(n; interval = SymmetricInterval())

            @test eltype(TanhSinhQuadrature(Float32, n; IT = Float32)) == Float32
            @test eltype(tanh_sinh_nodes(Float32, n; IT = Float32)) == Float32
            @test eltype(tanh_sinh_nodes(Float32, n; IT = Float32, interval = SymmetricInterval())) ==
                  Float32
            @test eltype(tanh_sinh_weights(Float32, n; IT = Float32)) == Float32
            @test eltype(tanh_sinh_weights(Float32, n; IT = Float32, interval = SymmetricInterval())) ==
                  Float32
        end
    end

    # Gauss-Legendre is the better rule for an integrand that is smooth up to and including
    # the endpoints; tanh-sinh wins decisively as soon as there is an endpoint singularity.
    @testset "$(rpad("comparison with Gauss-Legendre",60))" begin
        ts = TanhSinhQuadrature(3)                  # 51 nodes
        gl = GaussLegendreQuadrature(51; fast = true)

        @test abs(gl(x -> 1/(1+x^2)) - atan(1)) ≤ abs(ts(x -> 1/(1+x^2)) - atan(1))
        @test abs(ts(log) + 1) < abs(gl(log) + 1) / 1E6
    end
end
