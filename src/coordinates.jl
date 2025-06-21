"""
    transform(time, pos, in, out)

Transform coordinates from `in` coordinate system to `out` coordinate system.

Note: `pos` must be of shape (3,) for single point or (3, n) for multiple points

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
    size(pos, 1) != 3 && error("Position array must of shape (3, n), got size ", size(pos))

    # Prepare call arguments
    pos_in = arrf(pos)
    pos_out = similar(pos_in)
    ntime, iyear, idoy, ut = prepare_time(time)
    sys_in = coord_sys(in)
    sys_out = coord_sys(out)
    coord_trans_vec1!(ntime, sys_in, sys_out, iyear, idoy, ut, pos_in, pos_out)

    return isa(time, AbstractVector) ? pos_out : vec(pos_out)
end

transform(time, pos, inout) = transform(time, pos, parse_coord_transform(inout)...)

(::Type{S})(time, pos::CoordinateVector) where {S<:AbstractCoordinateSystem} =
    S(transform(time, pos, pos.sym, S))


# transform(time, pos::CoordinateVector, out) = transform(time, pos, pos.sym, out)
# transform(time, pos::CoordinateVector, ::Type{S}) where {S<:AbstractCoordinateSystem} =
#     S(transform(time, pos, pos.sym, S))
