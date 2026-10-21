"""
    coord_trans(in, out, time, pos)

Transform `pos` (shape `(3,)` or `(3, n)`, one time per point) between IRBEM coordinate systems named `in` and `out`.
"""
coord_trans(in, out, time, pos) = _coord_trans(coord_sys(in), coord_sys(out), time, pos)

function _coord_trans(in, out, time, pos)
    size(pos, 1) == 3 || throw(DimensionMismatch("positions must be of shape (3,) or (3, n), got size $(size(pos))"))
    iyear, idoy, ut = decompose_time(time)
    pos_in = convert(Array{Float64}, pos)
    n = size(pos_in, 2)
    length(ut) == n || throw(DimensionMismatch("got $(length(ut)) time(s) for $n position(s)"))
    pos_out = similar(pos_in)
    coord_trans_vec1!(Int32(n), in, out, iyear, idoy, ut, pos_in, pos_out)
    return pos_out
end
