"""
    transform(from => to, x, t; dim)
    transform(from => to, A; dim)
    transform(to, x::CoordinateVector, t)

Transform `x` from coordinate system `from` to `to` at time(s) `t` with IRBEM, with the call shapes and inputs of `GeoCotrans.transform`.
Coordinate systems and component orders: see README "Coordinates".
"""
function transform end

transform(p::Pair, x, t; dim = nothing) = _transform(_rawsys(p.first), _rawsys(p.second), x, t, dim)
transform(p::Pair, A::AbstractArray; dim = nothing) = transform(p, A, times(A, dim); dim)

function transform(to, x::CoordinateVector, t)
    F, R = _csys(to)
    y = _transform(_sys(x), _sys((F, R)), SVector{3, Float64}(x), t, nothing)
    return CoordinateVector{typeof(F), typeof(R)}(y...)
end

function transform(p::Pair, x::CoordinateVector{F}, t) where {F}
    _csys(p.first) in ((F(), Cartesian3()), getcsys(x)) || throw(ArgumentError("$(typeof(x)) is not in $(p.first)"))
    return transform(p.second, x, t)
end

function _transform((from, rin), (to, rout), x::AbstractVector{<:Real}, t, _)
    iyear, idoy, ut = decompose_time_s(_scalar(t))
    pos_out = MVector{3, Float64}(undef)
    coord_trans_vec1!(Int32(1), from, to, Ref(iyear), Ref(idoy), Ref(ut), _to_irbem(rin, SVector{3, Float64}(unwrap(x))), pos_out)
    y = _from_irbem(rout, SVector(pos_out))
    return x isa StaticArray ? y : rebuild(x, Vector(y))
end

function _transform((from, rin), (to, rout), x::AbstractVector{<:AbstractVector}, t, _)
    pos_in = Matrix{Float64}(undef, 3, length(x))
    for (j, v) in enumerate(unwrap(x))
        pos_in[:, j] = _to_irbem(rin, SVector{3, Float64}(v))
    end
    pos_out = _trans(from, to, t, pos_in)
    return rebuild(x, [_from_irbem(rout, SVector{3}(c)) for c in eachcol(pos_out)])
end

function _transform((from, rin), (to, rout), A::AbstractMatrix, t, dim)
    dim = _sampledim(A, dim)
    U = unwrap(A)
    size(U, 3 - dim) == 3 || throw(DimensionMismatch("expected 3 components along dimension $(3 - dim), got size $(size(U))"))
    # a 3×n Float64 matrix goes to Fortran as is
    pos_in = dim == 2 && _isidentity(rin) ? _arrf(U) : _stack3(rin, U, dim)
    pos_out = _trans(from, to, t, pos_in)
    if !_isidentity(rout)
        for c in eachcol(pos_out)
            c .= _from_irbem(rout, SVector{3}(c))
        end
    end
    return rebuild(A, dim == 2 ? pos_out : permutedims(pos_out))
end

function _stack3(rep, U, dim)
    out = Matrix{Float64}(undef, 3, size(U, dim))
    for j in axes(out, 2)
        x = dim == 2 ? SA[U[1, j], U[2, j], U[3, j]] : SA[U[j, 1], U[j, 2], U[j, 3]]
        out[:, j] = _to_irbem(rep, SVector{3, Float64}(x))
    end
    return out
end

function _trans(from, to, t, pos_in)
    n = size(pos_in, 2)
    iyear, idoy, ut = _decompose(t, n)
    pos_out = similar(pos_in)
    coord_trans_vec1!(Int32(n), from, to, iyear, idoy, ut, pos_in, pos_out)
    return pos_out
end

# `t` follows broadcast rules, as in GeoCotrans
function _decompose(t, n)
    bt = Base.broadcastable(t)
    isempty(axes(bt)) || return length(bt) == n ? decompose_time(bt) : throw(DimensionMismatch("got $(length(bt)) time(s) for $n sample(s)"))
    iyear, idoy, ut = decompose_time_s(bt[])
    return fill(iyear, n), fill(idoy, n), fill(ut, n)
end

function _scalar(t)
    bt = Base.broadcastable(t)
    isempty(axes(bt)) && return bt[]
    length(bt) == 1 || throw(DimensionMismatch("got $(length(bt)) time(s) for a single position"))
    return only(bt)
end

_sampledim(A, dim) = dimnum(A, dim)
_sampledim(A, ::Nothing) = something(tdimnum(A), 2)

_arrf(x) = convert(Array{Float64}, x)
