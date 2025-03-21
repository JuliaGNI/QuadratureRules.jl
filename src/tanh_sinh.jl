"""
tanh-sinh quadrature.
"""
function TanhSinhQuadrature(::Type{T}, n::Integer; IT=BigFloat) where {T}
    if n < 1
        throw(ErrorException("Tanh-Sinh quadrature is not defined for less than one level."))
    end

    h = IT(1)
    b = IT[IT(π)/2]
    c = IT[0]

    for l in 1:n
        k = 1
        h /= 2
        b ./= 2

        while true
            w = IT(π) / 2 * h * cosh(k * h) / ( cosh(IT(π) / 2 * sinh(k * h) ) )^2
            x = IT(1) / ( exp( IT(π) / 2 * sinh(k * h) ) * cosh( IT(π) / 2 * sinh(k * h) ) )

            T(w) == 0 && break
            T(abs(1-x)) == 1 && break
            T(abs(x-1)) == 1 && break

            push!(c, 1-x)
            push!(b, w)
            push!(c, x-1)
            push!(b, w)

            k += 1
            l == 1 || (k += 1)
        end
    end

    shift!(b, c)

    QuadratureRule(n, c, b, T)
end

TanhSinhQuadrature(n) = TanhSinhQuadrature(Float64, n)
