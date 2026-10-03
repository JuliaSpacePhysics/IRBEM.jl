"""
    transform(time, pos, in, out)
    transform(time, pos, in => out)
    transform(time, pos, "in2out")

Transform coordinates from `in` coordinate system to `out` coordinate system.

`pos` is of shape (3,) for a single point or (3, n) for `n` points, with one time per point.

# Example
```julia
using Dates
using IRBEM

time = DateTime(2020, 1, 1)
pos = [6.90274, -1.63624, 1.91669]

GSM(time, GEO(pos))
transform(time, pos, "GEO", "GSM")
transform(time, pos, "GEO" => "GSM")
transform(time, pos, "geo2gsm")
```
"""
function transform(time, pos, in, out)
    size(pos, 1) == 3 || throw(DimensionMismatch("positions must be of shape (3,) or (3, n), got size $(size(pos))"))
    iyear, idoy, ut = decompose_time(time)
    n = length(ut)
    npos = ndims(pos) == 1 ? 1 : size(pos, 2)
    n == npos || throw(DimensionMismatch("got $n time(s) but $npos position(s)"))
    pos_in = _arrf(pos)
    pos_out = similar(pos_in)
    coord_trans_vec1!(Int32(n), coord_sys(in), coord_sys(out), iyear, idoy, ut, pos_in, pos_out)
    return time isa AbstractVector ? pos_out : vec(pos_out)
end

transform(time, pos, inout) = transform(time, pos, parse_coord_transform(inout)...)

_arrf(x) = convert(Array{Float64}, x)
_arrf(x::StaticArray) = Float64.(SArray(x))

(::Type{S})(time, pos::CoordinateVector) where {S <: AbstractCoordinateSystem} =
    S(transform(time, pos, pos.sym, S))
