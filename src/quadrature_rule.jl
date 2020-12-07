
struct QuadratureRule{T,N}
    order::Int
    nodes::Vector{T}
    weights::Vector{T}
end

function QuadratureRule(order, nodes::AbstractVector{T}, weights::AbstractVector{T}, TT=T) where {T}
    @assert length(nodes) == length(weights)
    QuadratureRule{TT, length(nodes)}(order, nodes, weights)
end

"""
Integrate a function `f(x)` over the interval [0,1] using the quadrature.
"""
function (quad::QuadratureRule)(f::Function)
    mapreduce(i -> quad.weights[i] * f(quad.nodes[i]), +, eachindex(quad))
end


nnodes(::QuadratureRule{T,N}) where {T,N} = N
nodes(quad::QuadratureRule) = quad.nodes
order(quad::QuadratureRule) = quad.order
weights(quad::QuadratureRule) = quad.weights

Base.eachindex(quad::QuadratureRule) = eachindex(quad.nodes, quad.weights)

Base.hash(quad::QuadratureRule, h::UInt) = hash(quad.order, hash(quad.nodes, hash(quad.weights, hash(:QuadratureRule, h))))

Base.:(==)(quad1::QuadratureRule, quad2::QuadratureRule) = (
                   quad1.order == quad2.order
                && quad1.nodes == quad2.nodes
                && quad1.weights == quad2.weights)

Base.isequal(quad1::QuadratureRule{T1,N1}, quad2::QuadratureRule{T2,N2}) where {T1,N1,T2,N2} = (
                quad1 == quad2
                && T1 == T2
                && N1 == N2)

Base.isapprox(quad1::QuadratureRule, quad2::QuadratureRule; kwargs...) = (
                            quad1.order == quad2.order
                && isapprox(quad1.nodes,   quad2.nodes;   kwargs...)
                && isapprox(quad1.weights, quad2.weights; kwargs...))

Base.eltype(::QuadratureRule{T}) where {T} = T
