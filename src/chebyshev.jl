
"""
Chebyshev points of the first (`Val(1)`) or second (`Val(2)`) kind on the interval [-1,+1].
"""
function chebyshev_points(::Type{T}, s::Integer, ::Val{1}) where {T}
    T[ @big sin( π*(s-2i+1) / (2s) ) for i in s:-1:1 ]
end

function chebyshev_points(::Type{T}, s::Integer, ::Val{2}) where {T}
    if s == 1
        throw(ErrorException("Chebyshev points of the second kind are not defined for one point."))
    end

    T[ @big cos( π * (i-1) / (s-1) ) for i in s:-1:1 ]
end

chebyshev_points(s, kind) = chebyshev_points(Float64, s, Val(kind))

"""
Chebyshev nodes of the first (`Val(1)`) or second (`Val(2)`) kind on the interval [0,1].
"""
function chebyshev_nodes(::Type{T}, s::Integer, kind::Val) where {T}
    shift_nodes(chebyshev_points(T, s, kind))
end

chebyshev_nodes(s, kind) = chebyshev_nodes(Float64, s, Val(kind))


"Gauss-Chebyshev points on the interval [-1,+1]."
gauss_chebyshev_points(::Type{T}, s::Integer) where {T} = chebyshev_points(T, s, Val(1))
gauss_chebyshev_points(s) = gauss_chebyshev_points(Float64, s)

"Gauss-Chebyshev nodes on the interval [0,1]."
gauss_chebyshev_nodes(::Type{T}, s::Integer) where {T} = chebyshev_nodes(T, s, Val(1))
gauss_chebyshev_nodes(s) = gauss_chebyshev_nodes(Float64, s)

"Lobatto-Chebyshev points on the interval [-1,+1]."
lobatto_chebyshev_points(::Type{T}, s::Integer) where {T} = chebyshev_points(T, s, Val(2))
lobatto_chebyshev_points(s) = lobatto_chebyshev_points(Float64, s)

"Lobatto-Chebyshev nodes on the interval [0,1]."
lobatto_chebyshev_nodes(::Type{T}, s::Integer) where {T} = chebyshev_nodes(T, s, Val(2))
lobatto_chebyshev_nodes(s) = lobatto_chebyshev_nodes(Float64, s)

"""
Gauss-Chebyshev quadrature.
"""
function GaussChebyshevQuadrature(::Type{T}, s::Integer) where {T}
    local tj::BigFloat = 0
    local th::BigFloat = 0

    c = chebyshev_nodes(BigFloat, s, Val(1))
    b = zero(c)
    for i in eachindex(b)
        tj = 0
        th = @big π * (2i-1) / (2s)
        for j in 1:div(s,2)
            tj += @big cos(2j*th) / (4j^2 - 1)
        end
        b[i] = @big (1 - 2tj) / s
    end

    QuadratureRule(2s, c, b, T)
end

GaussChebyshevQuadrature(s) = GaussChebyshevQuadrature(Float64, s)


"""
Lobatto-Chebyshev quadrature.
"""
function LobattoChebyshevQuadrature(::Type{T}, s::Integer) where {T}
    if s == 1
        throw(ErrorException("Lobatto-Chebyshev quadrature is not defined for one stage."))
    end

    local tj::BigFloat = 0
    local th::BigFloat = 0

    c = chebyshev_nodes(BigFloat, s, Val(2))
    b = zero(c)
    for i in eachindex(b)
        tj = 0
        th = @big π * i / (s+1)
        for j in 1:div(s+1,2)
            tj += @big sin((2j-1) * th) / (2j - 1)
        end
        b[i] = @big (2tj * sin(th)) / (s+1)
    end

    QuadratureRule(2s, c, b, T)
end

LobattoChebyshevQuadrature(s) = LobattoChebyshevQuadrature(Float64, s)


ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{1}) where {T} = GaussChebyshevQuadrature(T,s)
ChebyshevQuadrature(::Type{T}, s::Integer, ::Val{2}) where {T} = LobattoChebyshevQuadrature(T,s)

ChebyshevQuadrature(s, kind) = ChebyshevQuadrature(Float64, s, Val(kind))
