
_big(x) = x
_big(x::Number) = big(x)
_big(x::String) = parse(BigFloat, x)

function _big(x::Expr)
    y = x
    y.args .= _big.(y.args)
    return  y
end

macro big(x)
    return esc(_big(x))
end


"Shift and scale nodes from the interval [-1,+1] to the interval [0,1]."
function shift_nodes(c)
    (c .+ 1) ./ 2
end

"Scale weights from the interval [-1,+1] to the interval [0,1]."
function scale_weights(b)
    b ./ 2
end

"Scale nodes and weights from the interval [-1,+1] to the interval [0,1]."
function shift!(b,c)
    b .= b ./ 2
    c .= (c .+ 1) ./ 2
end

"Scale nodes and weights from the interval [0,1] to the interval [-1,+1]."
function unshift!(b,c)
    b .= b .* 2
    c .= 2 .* c .- 1
end
