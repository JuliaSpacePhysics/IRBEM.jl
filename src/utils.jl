const BADDATA = -1.0e31 # IRBEM's output fill value
const MAGINPUT_MISSING = -9999.0
const MAGINPUT_FIELDS = fieldnames(MagInput) # rows 1-17 of IRBEM's maginput(25, ntime); rows 18-25 are reserved

_nan(x::Float64) = x == BADDATA ? NaN : x
_nan(x::Array{Float64}) = map!(_nan, x, x)
_nan(x::Ref) = _nan(x[])
_nan(x) = x

_squeeze(x::AbstractVector) = only(x)
_squeeze(x::AbstractArray) = dropdims(x; dims = ndims(x))

# A scalar time gives single-point outputs, whatever the number of points
_output(nt, single) = single ? map(_squeeze ∘ _nan, nt) : map(_nan, nt)

_out(p, dims...) = Array{Float64}(undef, dims..., p.ntime)

"""
    get_datetime(X::AbstractDict)

Extract datetime from input dictionary X.
Supports 'dateTime' or 'Time' keys with DateTime or String values.
"""
function get_datetime(X::AbstractDict)
    dt_val = get(X, "dateTime", get(X, "Time", nothing))
    return !isnothing(dt_val) ? Dates.DateTime.(dt_val) : error("No date/time information found in input dictionary. Expected 'dateTime' or 'Time' key.")
end

const CoordVectors = Union{CoordinateVector, AbstractVector{<:CoordinateVector}}
const CoordSpec = Union{
    AbstractString, Symbol, Integer, AbstractReferenceFrame, Type{<:AbstractReferenceFrame},
    typeof(GDZ), Tuple{AbstractReferenceFrame, AbstractRepresentation},
}

"""
    prepare_irbem(args...; kw...) -> (p, single)

Convert user inputs into the `(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput)`
arguments shared by IRBEM routines, as a `NamedTuple` `p` in that order.
`single` is `true` when `time` is a scalar.
"""
prepare_irbem(time, x, coord = "GDZ", maginput = (;); kext = KEXT[], options = OPTIONS[]) =
    _prepare(time, prepare_loc(x), _rawsys(coord), maginput, kext, options)

prepare_irbem(time, x::CoordVectors, maginput = (;); kext = KEXT[], options = OPTIONS[]) =
    _prepare(time, prepare_loc(x), _sys(x), maginput, kext, options)

function prepare_irbem(time, x::CoordVectors, coord::CoordSpec, maginput = (;); kw...)
    _csys(coord) == _getcsys(x) || throw(ArgumentError("coord = $coord contradicts the coordinate system of the input positions"))
    return prepare_irbem(time, x, maginput; kw...)
end

prepare_irbem(model::MagneticField, X::AbstractDict, maginput = (;)) =
    _prepare(get_datetime(X), prepare_loc(X["x1"], X["x2"], X["x3"]), (model.sysaxes, nothing), maginput, model.kext, model.options)

function _prepare(time, loc, (sysaxes, rep), maginput, kext, options)
    x1, x2, x3 = _to_irbem_loc(rep, loc)
    iyear, idoy, ut = decompose_time(time)
    n = length(ut)
    length(x1) == n || throw(DimensionMismatch("got $n time(s) but $(length(x1)) position(s)"))
    p = (;
        ntime = Int32(n), kext = parse_kext(kext), options = prepare_options(options), sysaxes = Int32(sysaxes),
        iyear, idoy, ut, x1, x2, x3, maginput = prepare_maginput(maginput, n),
    )
    return p, !(time isa AbstractVector)
end

# Inputs at point `i`, for IRBEM's single-point routines
function _point(p, i)
    return (
        p.kext, p.options, p.sysaxes, Ref(p.iyear, i), Ref(p.idoy, i), Ref(p.ut, i),
        Ref(p.x1, i), Ref(p.x2, i), Ref(p.x3, i), Ref(p.maginput, 25 * (i - 1) + 1),
    )
end

# Point `i` of an output whose last dimension indexes points
_at(x::Array, i) = Ref(x, stride(x, ndims(x)) * (i - 1) + 1)

# Call IRBEM's single-point routine `f!` at each point
function _each_point!(f!, p, extra, nt)
    for i in 1:p.ntime
        f!(_point(p, i)..., extra..., map(x -> _at(x, i), nt)...)
    end
    return nt
end

function _single((p, _))
    p.ntime == 1 || throw(ArgumentError("expected a single time and position, got $(p.ntime); broadcast over points instead"))
    return _point(p, 1)
end

function decompose_time_s(dt::DateTime)
    iyear = Int32(year(dt))
    idoy = Int32(dayofyear(dt))
    ut = Float64(hour(dt) * 3600 + minute(dt) * 60 + second(dt) + millisecond(dt) / 1000)
    return iyear, idoy, ut
end

decompose_time_s(dt) = decompose_time_s(DateTime(dt))

function decompose_time(ts::Union{AbstractVector, Tuple})
    n = length(ts)
    iyear = Vector{Int32}(undef, n)
    idoy = Vector{Int32}(undef, n)
    ut = Vector{Float64}(undef, n)
    for (i, t) in enumerate(ts)
        iyear[i], idoy[i], ut[i] = decompose_time_s(t)
    end
    return iyear, idoy, ut
end

decompose_time(t) = decompose_time((t,))

_vecf(x::Number) = [Float64(x)]
_vecf(x) = convert(Vector{Float64}, x)

prepare_loc(x1, x2, x3) = (_vecf(x1), _vecf(x2), _vecf(x3))
function prepare_loc(x)
    length(x) == 3 || throw(DimensionMismatch("a position must have 3 components, got $(length(x))"))
    return prepare_loc(x[1], x[2], x[3])
end
function prepare_loc(x::AbstractMatrix)
    size(x, 1) == 3 || throw(DimensionMismatch("positions must be a 3×n matrix, got size $(size(x))"))
    return prepare_loc(view(x, 1, :), view(x, 2, :), view(x, 3, :))
end
prepare_loc(x::AbstractVector{<:AbstractVector}) = prepare_loc(getindex.(x, 1), getindex.(x, 2), getindex.(x, 3))

function prepare_options(options)
    length(options) == 5 || throw(ArgumentError("options must have 5 elements, got $(length(options))"))
    # A fresh copy: some routines (e.g. landi2lstar) overwrite `options`
    return MVector{5, Int32}(options)
end

"""
    prepare_maginput(maginput, n)

Build IRBEM's 25×n `maginput` array from a `MagInput`, `NamedTuple` or `AbstractDict`.
Each value is a scalar (repeated for every point) or a vector of length `n`.
"""
function prepare_maginput(maginput, n)
    out = fill(MAGINPUT_MISSING, 25, n)
    for (k, v) in _pairs(maginput)
        i = findfirst(==(Symbol(k)), MAGINPUT_FIELDS)
        isnothing(i) && throw(ArgumentError("Unknown magnetic field input $k. Valid inputs are $(join(MAGINPUT_FIELDS, ", "))"))
        out[i, :] .= v
    end
    return out
end

_pairs(m::MagInput) = (f => getfield(m, f) for f in MAGINPUT_FIELDS)
_pairs(m) = pairs(m)

parse_kext(kext::Integer) = Int32(kext)
parse_kext(kext::ExternalFieldModel) = Int32(kext)
const _EXT_MODELS = uppercase.(string.(instances(ExternalFieldModel)))
function parse_kext(kext)
    idx = findfirst(==(uppercase(string(kext))), _EXT_MODELS)
    isnothing(idx) && throw(ArgumentError("Unknown external field model: $kext. Valid models are $(join(instances(ExternalFieldModel), ", "))"))
    return Int32(idx - 1)
end

# Position gives IRBEM's `sysaxes` code (0-based)
const SYSAXES = (:GDZ, :GEO, :GSM, :GSE, :SM, :GEI, :MAG, :SPH, :RLL, :HEE, :HAE, :HEEQ, :TOD, :J2000, :TEME)

"IRBEM `sysaxes` code of a coordinate system name (case-insensitive), code, or Cartesian frame."
function coord_sys(x::Union{Symbol, AbstractString})
    idx = findfirst(==(Symbol(x)), SYSAXES)
    isnothing(idx) && (idx = findfirst(==(Symbol(uppercase(String(x)))), SYSAXES))
    isnothing(idx) && throw(ArgumentError("Unknown coordinate system: $x. Choose from $(join(SYSAXES, ", "))."))
    return Int32(idx - 1)
end
coord_sys(x::Integer) = Int32(x)
coord_sys(::Type{F}) where {F <: AbstractReferenceFrame} = coord_sys(nameof(F))
coord_sys(f::AbstractReferenceFrame) = coord_sys(typeof(f))

# `(frame, representation)` of each `sysaxes` code; TOD is GEI in IRBEM, RLL has no GeoCotrans representation
const SYSAXES_CSYS = (
    GDZ(), (GEO(), Cartesian3()), (GSM(), Cartesian3()), (GSE(), Cartesian3()), (SM(), Cartesian3()),
    (GEI(), Cartesian3()), (MAG(), Cartesian3()), (GEO(), Spherical()), nothing,
    (HEE(), Cartesian3()), (HAE(), Cartesian3()), (HEEQ(), Cartesian3()), (GEI(), Cartesian3()),
    (J2000(), Cartesian3()), (TEME(), Cartesian3()),
)

function _csys(x::Union{Symbol, AbstractString, Integer})
    c = coord_sys(x)
    csys = 0 <= c < length(SYSAXES) ? SYSAXES_CSYS[c + 1] : nothing
    return @something csys throw(ArgumentError("$x has no typed coordinate equivalent"))
end
_csys(f::AbstractReferenceFrame) = (f, Cartesian3())
_csys(::Type{F}) where {F <: AbstractReferenceFrame} = (F(), Cartesian3())
_csys(::typeof(GDZ)) = GDZ()
_csys(c::Tuple{AbstractReferenceFrame, AbstractRepresentation}) = c
_csys(x) = throw(ArgumentError("Unknown coordinate system $(repr(x))"))

_getcsys(x::CoordinateVector) = getcsys(x)
_getcsys(::AbstractVector{<:CoordinateVector{F, R}}) where {F, R} = (F(), R())

# `(sysaxes, representation)` to call IRBEM with; names and codes pass components in IRBEM's order
_sys(x::Union{Symbol, AbstractString, Integer}) = (coord_sys(x), nothing)
_sys(x::CoordVectors) = _sys(_getcsys(x))
_sys(x) = _sys(_csys(x))
_sys((f, r)::Tuple{AbstractReferenceFrame, AbstractRepresentation}) = (coord_sys(f), r)
_sys((_, r)::Tuple{GEO, Geodetic}) = (Int32(0), r)
_sys((f, _)::Tuple{AbstractReferenceFrame, Geodetic}) = throw(ArgumentError("IRBEM has no geodetic $(nameof(typeof(f))) coordinates"))

# Raw components carry no representation, so `GDZ` and "GDZ" would read the same array in opposite orders:
# only names and codes (IRBEM's order) or frames (Cartesian) describe them; typed representations need a `CoordinateVector`.
_rawsys(x::Union{Symbol, AbstractString, Integer}) = _sys(x)
_rawsys(f::Union{AbstractReferenceFrame, Type{<:AbstractReferenceFrame}}) = (coord_sys(f), nothing)
_rawsys(x) = throw(ArgumentError("raw components in $(repr(x)) are ambiguous; wrap them in a `CoordinateVector` (e.g. `GDZ(lat, lon, alt)`) or use an IRBEM name such as \"GDZ\" (IRBEM's component order)"))

# Spherical positions go to IRBEM as Cartesian in their frame
_isidentity(rep) = rep isa Union{Nothing, Cartesian3}
_to_irbem(_, x) = x
_to_irbem(::Spherical, x) = sph2car(x...)
function _to_irbem(::Geodetic, (lat, lon, alt))
    abs(lat) > 90 && throw(ArgumentError("geodetic latitude $lat is outside [-90, 90]; `GDZ` takes (lat, lon, alt)"))
    return SA[alt, lat, lon]
end
_from_irbem(_, x) = x
_from_irbem(::Spherical, x) = car2sph(x...)
_from_irbem(::Geodetic, (alt, lat, lon)) = SA[lat, lon, alt]

_to_irbem_loc(rep, loc) = _isidentity(rep) ? loc : _unzip3(map((a, b, c) -> _to_irbem(rep, SA[a, b, c]), loc...))
_unzip3(v) = (getindex.(v, 1), getindex.(v, 2), getindex.(v, 3))

function clean_posit!(posit::Array{T, 3}, Blocal::Matrix{T}, Nposit) where {T}
    for (i, n) in enumerate(Nposit)
        posit[:, (n + 1):end, i] .= T(NaN)
        Blocal[(n + 1):end, i] .= T(NaN)
    end
    return
end
