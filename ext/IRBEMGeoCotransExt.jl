module IRBEMGeoCotransExt

using IRBEM: IRBEM, coord_sys, prepare_loc, _coord_trans
using GeoCotrans: CoordinateVector, Cartesian3, Spherical, Geodetic, GEO, sph2car
using StaticArrays: SVector, SA_F64

# IRBEM `sysaxes` code and components
_irbem(x::CoordinateVector{F, Cartesian3}) where {F} = coord_sys(String(nameof(F))), SVector{3, Float64}(x)
_irbem(x::CoordinateVector{F, Spherical}) where {F} = coord_sys(String(nameof(F))), SVector{3, Float64}(sph2car(x...))
_irbem(x::CoordinateVector{GEO, Geodetic}) = Int32(0), SA_F64[x[3], x[1], x[2]]

IRBEM._loc(x::CoordinateVector, coord) = IRBEM._loc([x], coord)
function IRBEM._loc(x::AbstractVector{<:CoordinateVector}, coord)
    isnothing(coord) || throw(ArgumentError("`coord` is given by the type of the positions"))
    ps = map(_irbem, x)
    sys = unique(first.(ps))
    length(sys) == 1 || throw(ArgumentError("positions must share one coordinate system, got $(length(sys))"))
    return prepare_loc(last.(ps)), only(sys)
end

function IRBEM._geo(x::CoordinateVector, t)
    sys, r = _irbem(x)
    return sys == 1 ? r : SVector{3, Float64}(_coord_trans(sys, Int32(1), t, r))
end

end
