import SymPyPythonCall
import SymPyPythonCall: simplify
import QuadratureRules: scale_weights, shift_nodes

"""
    symtype()

The symbolic element type, i.e. `Sym{T}` for `T` the Python type underlying `Sym(1)`.

Written as `typeof(Sym(1))` rather than spelled out, so that these tests do not depend on
which Python bridge SymPyPythonCall happens to be built on. This is the type with which
downstream packages such as RungeKutta.jl request symbolic coefficients.
"""
symtype() = typeof(SymPyPythonCall.Sym(1))

# Convert a symbolic value to a float by evaluating it. SymPyPythonCall's
# `convert(Float64, ::Sym)` uses `pyconvert`, which rejects unevaluated expressions such as
# `1/2 - sqrt(3)/6`, whereas the `N` path evaluates them. This is what makes the
# `symtype() ≈ Float64` comparisons below work, as their `isapprox` converts to `Float64`.
Base.convert(::Type{T}, x::SymPyPythonCall.Sym) where {T<:AbstractFloat} = T(SymPyPythonCall.N(x))

# Structural equality is not enough for algebraic numbers: `sqrt(3)^2 - 3` is zero but not
# syntactically so. Every exactness assertion below therefore simplifies the difference.
exactly(x, y) = iszero(simplify(x - y))
exactly(x::AbstractVector, y::AbstractVector) = length(x) == length(y) && all(exactly.(x, y))

"""
The checks that apply to every family: the four accessors return the requested element type
and the right number of values, they agree with the `Float64` result, the derived quantities
are consistent with the primary ones, and the weights sum to one *exactly*.

That last assertion is the one that distinguishes a genuinely exact result from a `BigFloat`
merely wrapped in a `Sym`, which is what these accessors used to return: a sum of rounded
weights is never symbolically equal to one.
"""
function test_symbolic_family(points, nodes, point_weights, weights, s)
    x = points(symtype(), s)
    c = nodes(symtype(), s)
    v = point_weights(symtype(), s)
    b = weights(symtype(), s)

    for y in (x, c, v, b)
        @test eltype(y) == symtype()
        @test length(y) == s
    end

    @test x ≈ points(Float64, s)
    @test c ≈ nodes(Float64, s)
    @test v ≈ point_weights(Float64, s)
    @test b ≈ weights(Float64, s)

    @test exactly(c, shift_nodes(x))
    @test exactly(b, scale_weights(v))

    @test exactly(sum(b), 1)
    @test exactly(sum(v), 2)
end

"""
The moment conditions ``\\sum_i b_i c_i^k = 1/(k+1)`` up to the degree of exactness `d`,
asserted as exact symbolic identities rather than to within a tolerance. This is the
strongest available statement that the symbolic path really computes the rule and does not
just carry floating point values around in symbolic clothing.
"""
function test_symbolic_moments(nodes, weights, s, d)
    c = nodes(symtype(), s)
    b = weights(symtype(), s)

    for k in 0:d
        @test exactly(sum(b .* c.^k), 1//(k+1))
    end
end

# Normalise the accessors that take an extra argument, so that they can be passed to the
# helpers above alongside the single-argument ones.
radau_points(endpoint)        = (T, s) -> radau_legendre_points(T, s, Val(endpoint))
radau_nodes(endpoint)         = (T, s) -> radau_legendre_nodes(T, s, Val(endpoint))
radau_point_weights(endpoint) = (T, s) -> radau_legendre_point_weights(T, s, Val(endpoint))
radau_weights(endpoint)       = (T, s) -> radau_legendre_weights(T, s, Val(endpoint))

chebyshev_pts(kind)          = (T, s) -> chebyshev_points(T, s, Val(kind))
chebyshev_nds(kind)          = (T, s) -> chebyshev_nodes(T, s, Val(kind))
chebyshev_pt_weights(kind)   = (T, s) -> chebyshev_point_weights(T, s, Val(kind))
chebyshev_wts(kind)          = (T, s) -> chebyshev_weights(T, s, Val(kind))


@testset "$(rpad("Symbolic evaluation",80))" begin

    # Every family other than tanh-sinh computes in the element type itself when that type
    # is not a floating point one, so that a computer algebra type gives exact results.
    # The number of stages is kept small throughout: SymPy resolves the nodes as radicals,
    # which becomes expensive beyond the quartic.

    @testset "$(rpad("Gauß-Legendre",60))" begin
        for s in 1:3
            test_symbolic_family(gauss_legendre_points, gauss_legendre_nodes,
                                 gauss_legendre_point_weights, gauss_legendre_weights, s)
            test_symbolic_moments(gauss_legendre_nodes, gauss_legendre_weights, s, 2s-1)
        end

        @test exactly(gauss_legendre_points(symtype(), 2), [-sqrt(SymPyPythonCall.Sym(3))/3,
                                                            sqrt(SymPyPythonCall.Sym(3))/3])
        @test exactly(gauss_legendre_nodes(symtype(), 2), [1//2 - sqrt(SymPyPythonCall.Sym(3))/6,
                                                           1//2 + sqrt(SymPyPythonCall.Sym(3))/6])
        @test exactly(gauss_legendre_weights(symtype(), 2), [1//2, 1//2])
    end

    @testset "$(rpad("Lobatto-Legendre",60))" begin
        for s in 2:4
            test_symbolic_family(lobatto_legendre_points, lobatto_legendre_nodes,
                                 lobatto_legendre_point_weights, lobatto_legendre_weights, s)
            test_symbolic_moments(lobatto_legendre_nodes, lobatto_legendre_weights, s, 2s-3)

            # the endpoints are included and exact, as they are for floating point types
            @test exactly(lobatto_legendre_points(symtype(), s)[begin], -1)
            @test exactly(lobatto_legendre_points(symtype(), s)[end], 1)
            @test exactly(lobatto_legendre_nodes(symtype(), s)[begin], 0)
            @test exactly(lobatto_legendre_nodes(symtype(), s)[end], 1)
        end

        @test exactly(lobatto_legendre_points(symtype(), 3), [-1, 0, 1])
        @test exactly(lobatto_legendre_weights(symtype(), 3), [1//6, 2//3, 1//6])
        @test exactly(lobatto_legendre_weights(symtype(), 4), [1//12, 5//12, 5//12, 1//12])

        @test_throws ErrorException lobatto_legendre_points(symtype(), 1)
        @test_throws ErrorException lobatto_legendre_nodes(symtype(), 1)
        @test_throws ErrorException lobatto_legendre_point_weights(symtype(), 1)
        @test_throws ErrorException lobatto_legendre_weights(symtype(), 1)
    end

    @testset "$(rpad("Radau-Legendre",60))" begin
        for endpoint in (:left, :right)
            for s in 1:3
                test_symbolic_family(radau_points(endpoint), radau_nodes(endpoint),
                                     radau_point_weights(endpoint), radau_weights(endpoint), s)
                test_symbolic_moments(radau_nodes(endpoint), radau_weights(endpoint), s, 2s-2)
            end
        end

        # the two variants are exact mirror images of each other
        for s in 1:3
            @test exactly(radau_legendre_points(symtype(), s, Val(:right)),
                          -reverse(radau_legendre_points(symtype(), s, Val(:left))))
            @test exactly(radau_legendre_weights(symtype(), s, Val(:right)),
                          reverse(radau_legendre_weights(symtype(), s, Val(:left))))
        end

        @test exactly(radau_legendre_nodes(symtype(), 2, Val(:left)), [0, 2//3])
        @test exactly(radau_legendre_weights(symtype(), 2, Val(:left)), [1//4, 3//4])
        @test exactly(radau_legendre_nodes(symtype(), 2, Val(:right)), [1//3, 1])
        @test exactly(radau_legendre_weights(symtype(), 2, Val(:right)), [3//4, 1//4])
    end

    @testset "$(rpad("Chebyshev",60))" begin
        # the closed forms are evaluated at SymPy's π rather than at a rounded one, so the
        # points come out as exact radicals
        for s in 1:5
            test_symbolic_family(chebyshev_pts(1), chebyshev_nds(1),
                                 chebyshev_pt_weights(1), chebyshev_wts(1), s)
            test_symbolic_moments(chebyshev_nds(1), chebyshev_wts(1), s, s-1)
        end

        for s in 2:5
            test_symbolic_family(chebyshev_pts(2), chebyshev_nds(2),
                                 chebyshev_pt_weights(2), chebyshev_wts(2), s)
            test_symbolic_moments(chebyshev_nds(2), chebyshev_wts(2), s, s-1)
        end

        @test exactly(chebyshev_points(symtype(), 2, Val(1)), [-sqrt(SymPyPythonCall.Sym(2))/2,
                                                                sqrt(SymPyPythonCall.Sym(2))/2])
        @test exactly(chebyshev_points(symtype(), 3, Val(2)), [-1, 0, 1])
        @test exactly(chebyshev_weights(symtype(), 3, Val(1)), [2//9, 5//9, 2//9])

        @test_throws ErrorException chebyshev_points(symtype(), 1, Val(2))
        @test_throws ErrorException chebyshev_nodes(symtype(), 1, Val(2))

        # the aliases dispatch to the same closed forms
        for s in 1:5
            @test gauss_chebyshev_points(symtype(), s) == chebyshev_points(symtype(), s, Val(1))
            @test gauss_chebyshev_weights(symtype(), s) == chebyshev_weights(symtype(), s, Val(1))
        end

        for s in 2:5
            @test lobatto_chebyshev_points(symtype(), s) == chebyshev_points(symtype(), s, Val(2))
            @test lobatto_chebyshev_weights(symtype(), s) == chebyshev_weights(symtype(), s, Val(2))
        end
    end

    @testset "$(rpad("Clenshaw-Curtis",60))" begin
        for s in 2:5
            test_symbolic_family(clenshaw_curtis_points, clenshaw_curtis_nodes,
                                 clenshaw_curtis_point_weights, clenshaw_curtis_weights, s)
            test_symbolic_moments(clenshaw_curtis_nodes, clenshaw_curtis_weights, s, s-1)
        end

        @test exactly(clenshaw_curtis_weights(symtype(), 3), [1//6, 2//3, 1//6])
        @test exactly(clenshaw_curtis_weights(symtype(), 5), [1//30, 4//15, 2//5, 4//15, 1//30])

        @test_throws ErrorException clenshaw_curtis_points(symtype(), 1)
        @test_throws ErrorException clenshaw_curtis_weights(symtype(), 1)
    end

    @testset "$(rpad("Quadrature rules",60))" begin
        # the constructors share the code paths of the accessors, so they are exact as well
        for quad in (GaussLegendreQuadrature(symtype(), 2),
                     LobattoLegendreQuadrature(symtype(), 3),
                     RadauLegendreQuadrature(symtype(), 2, Val(:left)),
                     RadauLegendreQuadrature(symtype(), 2, Val(:right)),
                     GaussChebyshevQuadrature(symtype(), 3),
                     LobattoChebyshevQuadrature(symtype(), 3),
                     ClenshawCurtisQuadrature(symtype(), 3))
            @test eltype(quad) == symtype()
            @test exactly(sum(weights(quad)), 1)
        end

        @test exactly(nodes(GaussLegendreQuadrature(symtype(), 2)),
                      gauss_legendre_nodes(symtype(), 2))
        @test exactly(weights(GaussLegendreQuadrature(symtype(), 2)),
                      gauss_legendre_weights(symtype(), 2))

        @test GaussLegendreQuadrature(symtype(), 2) ≈ GaussLegendreQuadrature(Float64, 2)
        @test LobattoLegendreQuadrature(symtype(), 3) ≈ LobattoLegendreQuadrature(Float64, 3)
    end

    @testset "$(rpad("Working arithmetic",60))" begin
        # `IT` still selects the arithmetic explicitly. Setting it to a floating point type
        # is the escape hatch for element types that merely wrap numbers and cannot be
        # evaluated exactly: the rule is then computed in `IT` and converted at the end.
        c = gauss_legendre_nodes(symtype(), 3; IT=BigFloat)

        @test eltype(c) == symtype()
        @test c ≈ gauss_legendre_nodes(Float64, 3)

        # what comes back that way are rounded numbers wearing symbolic clothing, so they
        # are not symbolically equal to the radicals the default arithmetic produces
        @test !exactly(c, gauss_legendre_nodes(symtype(), 3))

        # floating point element types are unaffected by the symbolic path
        @test gauss_legendre_nodes(Float64, 3) == gauss_legendre_nodes(Float64, 3; IT=BigFloat)
        @test gauss_legendre_nodes(BigFloat, 3) == gauss_legendre_nodes(BigFloat, 3; IT=BigFloat)
    end

    @testset "$(rpad("Tanh-Sinh",60))" begin
        # Tanh-Sinh is the one family without an exact variant: the number of nodes is not
        # given in advance but follows from where they stop being resolvable in `T`, which
        # a type that does not round cannot answer.
        @test_throws ArgumentError tanh_sinh_points(symtype(), 2)
        @test_throws ArgumentError tanh_sinh_nodes(symtype(), 2)
        @test_throws ArgumentError tanh_sinh_point_weights(symtype(), 2)
        @test_throws ArgumentError tanh_sinh_weights(symtype(), 2)
        @test_throws ArgumentError TanhSinhQuadrature(symtype(), 2)
    end

end
