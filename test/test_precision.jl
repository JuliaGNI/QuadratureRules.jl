
@testset "$(rpad("Working Precision",80))" begin

    # worst violation of the moment conditions, evaluated in BigFloat from the
    # nodes and weights the rule actually returned
    function moment_error(quad)
        c, b = big.(nodes(quad)), big.(weights(quad))
        maximum(abs(sum(b .* c.^k) - 1 / big(k+1)) for k in 0:order(quad)-1)
    end

    # the Radau family takes the endpoint as a third positional argument, so it is
    # bound here to give it the same signature as the other families
    radau_left_rule(T, s; kwargs...)  = RadauLegendreQuadrature(T, s, Val(:left); kwargs...)
    radau_right_rule(T, s; kwargs...) = RadauLegendreQuadrature(T, s, Val(:right); kwargs...)

    radau_left_points(T, s; kwargs...)  = radau_legendre_points(T, s, Val(:left); kwargs...)
    radau_right_points(T, s; kwargs...) = radau_legendre_points(T, s, Val(:right); kwargs...)

    radau_left_nodes(T, s; kwargs...)  = radau_legendre_nodes(T, s, Val(:left); kwargs...)
    radau_right_nodes(T, s; kwargs...) = radau_legendre_nodes(T, s, Val(:right); kwargs...)

    radau_left_weights(T, s; kwargs...)  = radau_legendre_weights(T, s, Val(:left); kwargs...)
    radau_right_weights(T, s; kwargs...) = radau_legendre_weights(T, s, Val(:right); kwargs...)

    radau_left_point_weights(T, s; kwargs...)  = radau_legendre_point_weights(T, s, Val(:left); kwargs...)
    radau_right_point_weights(T, s; kwargs...) = radau_legendre_point_weights(T, s, Val(:right); kwargs...)

    RULES = (GaussLegendreQuadrature, LobattoLegendreQuadrature, ClenshawCurtisQuadrature,
             GaussChebyshevQuadrature, LobattoChebyshevQuadrature,
             radau_left_rule, radau_right_rule)

    POINTS = (gauss_legendre_points, lobatto_legendre_points, gauss_chebyshev_points,
              lobatto_chebyshev_points, clenshaw_curtis_points,
              radau_left_points, radau_right_points)

    NODES = (gauss_legendre_nodes, lobatto_legendre_nodes, gauss_chebyshev_nodes,
             lobatto_chebyshev_nodes, clenshaw_curtis_nodes,
             radau_left_nodes, radau_right_nodes)

    WEIGHTS = (gauss_legendre_weights, lobatto_legendre_weights, gauss_chebyshev_weights,
               lobatto_chebyshev_weights, clenshaw_curtis_weights,
               radau_left_weights, radau_right_weights)

    POINT_WEIGHTS = (gauss_legendre_point_weights, lobatto_legendre_point_weights,
                     gauss_chebyshev_point_weights, lobatto_chebyshev_point_weights,
                     clenshaw_curtis_point_weights,
                     radau_left_point_weights, radau_right_point_weights)

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

        for weightfunction in WEIGHTS, s in 2:8
            b = weightfunction(T, s; IT=IT)
            @test eltype(b) == T
            @test length(b) == s
            @test all(bᵢ > 0 for bᵢ in b)
            @test sum(b) ≈ one(T)  atol = tol * eps(T)
        end

        for weightfunction in POINT_WEIGHTS, s in 2:8
            b = weightfunction(T, s; IT=IT)
            @test eltype(b) == T
            @test length(b) == s
            @test all(bᵢ > 0 for bᵢ in b)
            @test sum(b) ≈ 2one(T)  atol = 2 * tol * eps(T)
        end

        for kind in (1, 2), s in 2:8
            @test eltype(chebyshev_points(T, s, Val(kind); IT=IT)) == T
            @test eltype(chebyshev_nodes(T, s, Val(kind); IT=IT)) == T
        end

        for kind in (1, 2), s in 2:8
            @test eltype(ChebyshevQuadrature(T, s, Val(kind); IT=IT)) == T
        end

        # the Legendre rules additionally offer a double precision shortcut
        for quadrature in (GaussLegendreQuadrature, LobattoLegendreQuadrature,
                           radau_left_rule, radau_right_rule), s in 2:8
            @test eltype(quadrature(T, s; fast=true)) == T
        end

        # Tanh-sinh cannot join the loops above: it takes a level rather than a node count,
        # so neither its number of nodes nor `moment_error` — which needs a positive order —
        # applies. Its nodes must nevertheless meet the same structural requirements, and,
        # being the rule for endpoint singularities, must stay strictly inside the interval.
        for n in 1:5
            quad = TanhSinhQuadrature(T, n; IT=IT)
            c, b = nodes(quad), weights(quad)

            @test eltype(quad) == T
            @test nnodes(quad) == length(c) == length(b)
            @test all(0 < cᵢ < 1 for cᵢ in c)
            @test all(bᵢ > 0 for bᵢ in b)
            @test issorted(c, lt = <=)
            @test allunique(c)

            for accessor in (tanh_sinh_points, tanh_sinh_nodes)
                x = accessor(T, n; IT=IT)
                @test eltype(x) == T
                @test length(x) == nnodes(quad)
                @test issorted(x, lt = <=)
                @test allunique(x)
            end

            # its weights truncate with the nodes, so they too must match in length; they
            # are only asserted positive, the sum being subject to the truncation error
            for accessor in (tanh_sinh_point_weights, tanh_sinh_weights)
                w = accessor(T, n; IT=IT)
                @test eltype(w) == T
                @test length(w) == nnodes(quad)
                @test all(wᵢ > 0 for wᵢ in w)
                @test w == reverse(w)
            end
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
