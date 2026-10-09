# The helpers below carry names general enough to collide with those of another test file, so
# they live in a module of their own rather than in `Main`.
module TestSymbolic

using QuadratureRules
using Test

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

# Evaluate symbolic values numerically, so that they can be compared with the floating
# point results. `convert(Float64, ::Sym)` is of no use here: it goes through `pyconvert`,
# which rejects unevaluated expressions such as `1/2 - sqrt(3)/6`, whereas `N` evaluates
# them first. Converting explicitly at the comparison sites keeps this out of `Base.convert`,
# which neither this package nor its tests own.
tofloat(x::AbstractVector) = Float64[SymPyPythonCall.N(xᵢ) for xᵢ in x]

# Structural equality is not enough for algebraic numbers: `sqrt(3)^2 - 3` is zero but not
# syntactically so. Every exactness assertion below therefore simplifies the difference.
exactly(x, y) = iszero(simplify(x - y))
function exactly(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) && all(exactly.(x, y))
end

"""
The checks that apply to every family: the node and weight functions return the requested
element type and the right number of values on either interval, they agree with the `Float64`
result, the two intervals are consistent with each other, and the weights sum to one *exactly*.

That last assertion is the one that distinguishes a genuinely exact result from a `BigFloat`
merely wrapped in a `Sym`, which is what these accessors used to return: a sum of rounded
weights is never symbolically equal to one.
"""
function test_symbolic_family(nodefunction, weightfunction, s)
    x = nodefunction(symtype(), s; interval = SymmetricInterval())
    c = nodefunction(symtype(), s)
    v = weightfunction(symtype(), s; interval = SymmetricInterval())
    b = weightfunction(symtype(), s)

    for y in (x, c, v, b)
        @test eltype(y) == symtype()
        @test length(y) == s
    end

    @test tofloat(x) ≈ nodefunction(Float64, s; interval = SymmetricInterval())
    @test tofloat(c) ≈ nodefunction(Float64, s)
    @test tofloat(v) ≈ weightfunction(Float64, s; interval = SymmetricInterval())
    @test tofloat(b) ≈ weightfunction(Float64, s)

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
function test_symbolic_moments(nodefunction, weightfunction, s, d)
    c = nodefunction(symtype(), s)
    b = weightfunction(symtype(), s)

    for k in 0:d
        @test exactly(sum(b .* c .^ k), 1//(k+1))
    end
end

# Normalise the accessors that take an extra positional argument, so that they can be passed
# to the helpers above alongside the single-argument ones. The `interval` keyword is forwarded.
function radau_nds(endpoint)
    (T, s; kwargs...) -> radau_legendre_nodes(T, s, Val(endpoint); kwargs...)
end
function radau_wts(endpoint)
    (T, s; kwargs...) -> radau_legendre_weights(T, s, Val(endpoint); kwargs...)
end

chebyshev_nds(kind) = (T, s; kwargs...) -> chebyshev_nodes(T, s, Val(kind); kwargs...)
chebyshev_wts(kind) = (T, s; kwargs...) -> chebyshev_weights(T, s, Val(kind); kwargs...)

@testset "$(rpad("Symbolic evaluation",80))" begin

    # Every family other than tanh-sinh computes in the element type itself when that type
    # is not a floating point one, so that a computer algebra type gives exact results.
    # The number of stages is kept small throughout: SymPy resolves the nodes as radicals,
    # which becomes expensive beyond the quartic.

    @testset "$(rpad("Gauß-Legendre",60))" begin
        for s in 1:3
            test_symbolic_family(gauss_legendre_nodes, gauss_legendre_weights, s)
            test_symbolic_moments(gauss_legendre_nodes, gauss_legendre_weights, s, 2s-1)
        end

        @test exactly(gauss_legendre_nodes(symtype(), 2; interval = SymmetricInterval()),
            [-sqrt(SymPyPythonCall.Sym(3))/3,
                sqrt(SymPyPythonCall.Sym(3))/3])
        @test exactly(gauss_legendre_nodes(symtype(), 2),
            [1//2 - sqrt(SymPyPythonCall.Sym(3))/6,
                1//2 + sqrt(SymPyPythonCall.Sym(3))/6])
        @test exactly(gauss_legendre_weights(symtype(), 2), [1//2, 1//2])
    end

    @testset "$(rpad("Lobatto-Legendre",60))" begin
        for s in 2:4
            test_symbolic_family(lobatto_legendre_nodes, lobatto_legendre_weights, s)
            test_symbolic_moments(lobatto_legendre_nodes, lobatto_legendre_weights, s, 2s-3)

            # the endpoints are included and exact, as they are for floating point types
            @test exactly(lobatto_legendre_nodes(symtype(), s; interval = SymmetricInterval())[begin], -1)
            @test exactly(lobatto_legendre_nodes(symtype(), s; interval = SymmetricInterval())[end], 1)
            @test exactly(lobatto_legendre_nodes(symtype(), s)[begin], 0)
            @test exactly(lobatto_legendre_nodes(symtype(), s)[end], 1)
        end

        @test exactly(lobatto_legendre_nodes(symtype(), 3; interval = SymmetricInterval()), [
            -1, 0, 1])
        @test exactly(lobatto_legendre_weights(symtype(), 3), [1//6, 2//3, 1//6])
        @test exactly(lobatto_legendre_weights(symtype(), 4), [1//12, 5//12, 5//12, 1//12])

        @test_throws ErrorException lobatto_legendre_nodes(symtype(), 1; interval = SymmetricInterval())
        @test_throws ErrorException lobatto_legendre_nodes(symtype(), 1)
        @test_throws ErrorException lobatto_legendre_weights(symtype(), 1; interval = SymmetricInterval())
        @test_throws ErrorException lobatto_legendre_weights(symtype(), 1)
    end

    @testset "$(rpad("Radau-Legendre",60))" begin
        for endpoint in (:left, :right)
            for s in 1:3
                test_symbolic_family(radau_nds(endpoint), radau_wts(endpoint), s)
                test_symbolic_moments(radau_nds(endpoint), radau_wts(endpoint), s, 2s-2)
            end
        end

        # the two variants are exact mirror images of each other
        for s in 1:3
            @test exactly(
                radau_legendre_nodes(symtype(), s, Val(:right); interval = SymmetricInterval()),
                -reverse(radau_legendre_nodes(symtype(), s, Val(:left); interval = SymmetricInterval())))
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
        # nodes come out as exact radicals
        for s in 1:5
            test_symbolic_family(chebyshev_nds(1), chebyshev_wts(1), s)
            test_symbolic_moments(chebyshev_nds(1), chebyshev_wts(1), s, s-1)
        end

        for s in 2:5
            test_symbolic_family(chebyshev_nds(2), chebyshev_wts(2), s)
            test_symbolic_moments(chebyshev_nds(2), chebyshev_wts(2), s, s-1)
        end

        @test exactly(
            chebyshev_nodes(symtype(), 2, Val(1); interval = SymmetricInterval()),
            [-sqrt(SymPyPythonCall.Sym(2))/2,
                sqrt(SymPyPythonCall.Sym(2))/2])
        @test exactly(chebyshev_nodes(symtype(), 3, Val(2); interval = SymmetricInterval()), [
            -1, 0, 1])
        @test exactly(chebyshev_weights(symtype(), 3, Val(1)), [2//9, 5//9, 2//9])

        @test_throws ErrorException chebyshev_nodes(symtype(), 1, Val(2); interval = SymmetricInterval())
        @test_throws ErrorException chebyshev_nodes(symtype(), 1, Val(2))

        # the aliases dispatch to the same closed forms
        for s in 1:5
            @test gauss_chebyshev_nodes(symtype(), s; interval = SymmetricInterval()) ==
                  chebyshev_nodes(symtype(), s, Val(1); interval = SymmetricInterval())
            @test gauss_chebyshev_weights(symtype(), s) ==
                  chebyshev_weights(symtype(), s, Val(1))
        end

        for s in 2:5
            @test lobatto_chebyshev_nodes(symtype(), s; interval = SymmetricInterval()) ==
                  chebyshev_nodes(symtype(), s, Val(2); interval = SymmetricInterval())
            @test lobatto_chebyshev_weights(symtype(), s) ==
                  chebyshev_weights(symtype(), s, Val(2))
        end
    end

    @testset "$(rpad("Clenshaw-Curtis",60))" begin
        for s in 2:5
            test_symbolic_family(clenshaw_curtis_nodes, clenshaw_curtis_weights, s)
            test_symbolic_moments(clenshaw_curtis_nodes, clenshaw_curtis_weights, s, s-1)
        end

        @test exactly(clenshaw_curtis_weights(symtype(), 3), [1//6, 2//3, 1//6])
        @test exactly(clenshaw_curtis_weights(symtype(), 5), [
            1//30, 4//15, 2//5, 4//15, 1//30])

        @test_throws ErrorException clenshaw_curtis_nodes(symtype(), 1; interval = SymmetricInterval())
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
            ClenshawCurtisQuadrature(symtype(), 3),
            ChebyshevQuadrature(symtype(), 3, Val(1)),
            ChebyshevQuadrature(symtype(), 3, Val(2)))
            @test eltype(quad) == symtype()
            @test exactly(sum(weights(quad)), 1)
        end

        # the umbrella constructor forwards to the two Chebyshev rules unchanged
        @test ChebyshevQuadrature(symtype(), 3, Val(1)) ==
              GaussChebyshevQuadrature(symtype(), 3)
        @test ChebyshevQuadrature(symtype(), 3, Val(2)) ==
              LobattoChebyshevQuadrature(symtype(), 3)

        @test exactly(nodes(GaussLegendreQuadrature(symtype(), 2)),
            gauss_legendre_nodes(symtype(), 2))
        @test exactly(weights(GaussLegendreQuadrature(symtype(), 2)),
            gauss_legendre_weights(symtype(), 2))

        # evaluated numerically, a symbolic rule is the floating point rule
        evaluated(q) = QuadratureRule(order(q), tofloat(nodes(q)), tofloat(weights(q)))

        @test evaluated(GaussLegendreQuadrature(symtype(), 2)) ≈
              GaussLegendreQuadrature(Float64, 2)
        @test evaluated(LobattoLegendreQuadrature(symtype(), 3)) ≈
              LobattoLegendreQuadrature(Float64, 3)
    end

    @testset "$(rpad("Working arithmetic",60))" begin
        # `IT` still selects the arithmetic explicitly. Setting it to a floating point type
        # is the escape hatch for element types that merely wrap numbers and cannot be
        # evaluated exactly: the rule is then computed in `IT` and converted at the end.
        c = gauss_legendre_nodes(symtype(), 3; IT = BigFloat)

        @test eltype(c) == symtype()
        @test tofloat(c) ≈ gauss_legendre_nodes(Float64, 3)

        # what comes back that way are rounded numbers wearing symbolic clothing, so they
        # are not symbolically equal to the radicals the default arithmetic produces
        @test !exactly(c, gauss_legendre_nodes(symtype(), 3))

        # floating point element types are unaffected by the symbolic path
        @test gauss_legendre_nodes(Float64, 3) ==
              gauss_legendre_nodes(Float64, 3; IT = BigFloat)
        @test gauss_legendre_nodes(BigFloat, 3) ==
              gauss_legendre_nodes(BigFloat, 3; IT = BigFloat)
    end

    @testset "$(rpad("Tanh-Sinh",60))" begin
        # Tanh-Sinh is the one family without an exact variant: the number of nodes is not
        # given in advance but follows from where they stop being resolvable in `T`, which
        # a type that does not round cannot answer.
        @test_throws ArgumentError tanh_sinh_nodes(symtype(), 2; interval = SymmetricInterval())
        @test_throws ArgumentError tanh_sinh_nodes(symtype(), 2)
        @test_throws ArgumentError tanh_sinh_weights(symtype(), 2; interval = SymmetricInterval())
        @test_throws ArgumentError tanh_sinh_weights(symtype(), 2)
        @test_throws ArgumentError TanhSinhQuadrature(symtype(), 2)
    end
end

end # module TestSymbolic
