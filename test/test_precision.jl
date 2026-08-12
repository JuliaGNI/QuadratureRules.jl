
@testset "$(rpad("Working Precision",80))" begin

    # worst violation of the moment conditions, evaluated in BigFloat from the
    # nodes and weights the rule actually returned
    function moment_error(quad)
        c, b = big.(nodes(quad)), big.(weights(quad))
        maximum(abs(sum(b .* c.^k) - 1 / big(k+1)) for k in 0:order(quad)-1)
    end

    RULES = (GaussLegendreQuadrature, LobattoLegendreQuadrature, ClenshawCurtisQuadrature,
             GaussChebyshevQuadrature, LobattoChebyshevQuadrature)

    POINTS = (gauss_legendre_points, lobatto_legendre_points, gauss_chebyshev_points,
              lobatto_chebyshev_points, clenshaw_curtis_points)

    NODES = (gauss_legendre_nodes, lobatto_legendre_nodes, gauss_chebyshev_nodes,
             lobatto_chebyshev_nodes, clenshaw_curtis_nodes)

    # (T, IT, tolerance in multiples of eps(T)). Computing in a working precision
    # wider than T yields correctly rounded results, so the moment conditions hold
    # to better than eps(T). Setting IT = T is markedly worse: round-off in the
    # intermediate terms, and in the Newton refinement of the Legendre roots,
    # costs two to three orders of magnitude. This is why BigFloat is the default.
    COMBINATIONS = ((Float32, Float32, 1000),
                    (Float32, Float64, 1),
                    (Float32, BigFloat, 1),
                    (Float64, Float64, 1000),
                    (Float64, BigFloat, 1))

    for (T, IT, tol) in COMBINATIONS

        for quadrature in RULES, s in 2:8
            quad = quadrature(T, s; IT=IT)

            @test typeof(quad) <: QuadratureRule
            @test eltype(quad) == T
            @test nnodes(quad) == s
            @test length(nodes(quad)) == s
            @test length(weights(quad)) == s

            @test all(0 ≤ cᵢ ≤ 1 for cᵢ in nodes(quad))
            @test issorted(nodes(quad), lt = <)
            @test all(bᵢ > 0 for bᵢ in weights(quad))
            # the weight sum is the k = 0 moment, so it degrades with IT in the
            # same way and carries the same tolerance
            @test sum(weights(quad)) ≈ one(T)  atol = tol * eps(T)

            @test moment_error(quad) ≤ tol * eps(T)
        end

        for points in POINTS, s in 2:8
            x = points(T, s; IT=IT)
            @test eltype(x) == T
            @test length(x) == s
            @test all(-1 ≤ xᵢ ≤ 1 for xᵢ in x)
            @test issorted(x, lt = <)
        end

        for nodefunction in NODES, s in 2:8
            c = nodefunction(T, s; IT=IT)
            @test eltype(c) == T
            @test length(c) == s
            @test all(0 ≤ cᵢ ≤ 1 for cᵢ in c)
            @test issorted(c, lt = <)
        end

        for kind in (1, 2), s in 2:8
            @test eltype(chebyshev_points(T, s, Val(kind); IT=IT)) == T
            @test eltype(chebyshev_nodes(T, s, Val(kind); IT=IT)) == T
        end

        for kind in (1, 2), s in 2:8
            @test eltype(ChebyshevQuadrature(T, s, Val(kind); IT=IT)) == T
        end

        # the Legendre rules additionally offer a double precision shortcut
        for quadrature in (GaussLegendreQuadrature, LobattoLegendreQuadrature), s in 2:8
            @test eltype(quadrature(T, s; fast=true)) == T
        end
    end

    # the tabulated rules hold exact values and take no working precision
    for quadrature in (RiemannQuadratureLeft, RiemannQuadratureRight,
                       MidpointQuadrature, TrapezoidalQuadrature), T in (Float32, Float64, BigFloat)
        @test eltype(quadrature(T)) == T
        @test moment_error(quadrature(T)) == 0
    end

    # computing in a wider precision than requested must be at least as accurate
    # as computing in the target precision itself
    for quadrature in RULES, s in 2:8
        @test moment_error(quadrature(Float64, s; IT=BigFloat)) ≤
              moment_error(quadrature(Float64, s; IT=Float64)) + eps(Float64)
    end

end
