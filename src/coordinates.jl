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
transform(time, pos, "GEO", "GSM")
transform(time, pos, "GEO" => "GSM")
transform(time, pos, "geo2gsm")
```
"""
function transform(time, pos, in, out)
    size(pos, 1) != 3 && error("Position array must of shape (3, n), got size ", size(pos))

    # Prepare call arguments
    pos_in = Array{Float64}(pos)
    pos_out = similar(pos_in)
    ntime, iyear, idoy, ut = prepare_time(time)
    sys_in = coord_sys(in)
    sys_out = coord_sys(out)
    @ccall libirbem.coord_trans_vec1_(
        ntime::Ref{Int32}, sys_in::Ref{Int32}, sys_out::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        pos_in::Ptr{Float64}, pos_out::Ptr{Float64}
    )::Cvoid

    return isa(time, Array) ? pos_out : vec(pos_out)
end

transform(time, pos, pair) = transform(time, pos, pair[1], pair[2])
transform(time, pos, s::String) = transform(time, pos, parse_coord_transform(s)...)
