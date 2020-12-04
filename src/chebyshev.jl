
function GaussChebyshevQuadrature(::Type{T}, s::Integer) where {T}
    local tj::BigFloat = 0
    local th::BigFloat = 0

    c = shift_nodes(reverse([ @big sin( π*(s-2i+1) / (2s) ) for i in 1:s ], 1))
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


function LobattoChebyshevQuadrature(::Type{T}, s::Integer) where {T}
    local tj::BigFloat = 0
    local th::BigFloat = 0

    c = shift_nodes(reverse([ @big cos( π * (i-1) / (s-1) ) for i in 1:s ]))
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

ChebyshevQuadrature(s, ::Val{kind}) where {kind} = ChebyshevQuadrature(Float64, s, Val(kind))
