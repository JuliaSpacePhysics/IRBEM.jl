const BADDATA = -1.0e31 # IRBEM's output fill value
const MAGINPUT_MISSING = -9999.0
const MAGINPUT_FIELDS = fieldnames(MagInput) # rows 1-17 of IRBEM's maginput(25, ntime); rows 18-25 are reserved

_nan(x::Float64) = x == BADDATA ? NaN : x
_nan(x::Array{Float64}) = map!(_nan, x, x)
_nan(x::Ref) = _nan(x[])
_nan(x) = x

# Single-point inputs give single-point outputs
_squeeze(x::AbstractVector) = only(x)
_squeeze(x::AbstractArray) = dropdims(x; dims = ndims(x))

_output(nt, n) = (map(_nan, nt); n == 1 ? map(_squeeze, nt) : nt)

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

"""
    prepare_irbem(args...; kw...)

Convert user inputs into the `(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput)`
arguments shared by IRBEM routines, as a `NamedTuple` in that order.
"""
prepare_irbem(time, x, coord = "GDZ", maginput = (;); kext = KEXT[], options = OPTIONS[]) =
    _prepare(time, prepare_loc(x), coord_sys(coord), maginput, kext, options)

prepare_irbem(time, x::CoordVectors, maginput = (;); kext = KEXT[], options = OPTIONS[]) =
    _prepare(time, prepare_loc(x), coord_sys(x), maginput, kext, options)

prepare_irbem(model::MagneticField, X::AbstractDict, maginput = (;)) =
    _prepare(get_datetime(X), prepare_loc(X["x1"], X["x2"], X["x3"]), model.sysaxes, maginput, model.kext, model.options)

function _prepare(time, (x1, x2, x3), sysaxes, maginput, kext, options)
    iyear, idoy, ut = decompose_time(time)
    n = length(ut)
    length(x1) == n || throw(DimensionMismatch("got $n time(s) but $(length(x1)) position(s)"))
    return (;
        ntime = Int32(n), kext = parse_kext(kext), options = prepare_options(options), sysaxes = Int32(sysaxes),
        iyear, idoy, ut, x1, x2, x3, maginput = prepare_maginput(maginput, n),
    )
end

# Inputs at point `i`, for IRBEM's single-point routines
function _point(p, i)
    return (
        p.kext, p.options, p.sysaxes, Ref(p.iyear, i), Ref(p.idoy, i), Ref(p.ut, i),
        Ref(p.x1, i), Ref(p.x2, i), Ref(p.x3, i), Ref(p.maginput, 25 * (i - 1) + 1),
    )
end

function _single(p)
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

"IRBEM `sysaxes` code of coordinate system `x` (name, `Symbol`, type, instance or `CoordinateVector`)."
function coord_sys(x::Union{Symbol, AbstractString})
    idx = findfirst(==(Symbol(x)), COORD_SYSTEMS)
    isnothing(idx) && (idx = findfirst(==(Symbol(uppercase(String(x)))), COORD_SYSTEMS))
    isnothing(idx) && throw(ArgumentError("Unknown coordinate system: $x. Choose from $(join(COORD_SYSTEMS, ", "))."))
    return Int32(idx - 1)
end
coord_sys(x::Integer) = Int32(x)
coord_sys(::Type{S}) where {S <: AbstractCoordinateSystem} = coord_sys(nameof(S))
coord_sys(x::AbstractCoordinateSystem) = coord_sys(typeof(x))
coord_sys(v::CoordinateVector) = coord_sys(v.sym)
coord_sys(::AbstractVector{<:CoordinateVector{<:Any, C}}) where {C} = coord_sys(C)

parse_coord_transform(pair::Pair) = pair.first, pair.second
parse_coord_transform(s::Symbol) = parse_coord_transform(String(s))
# "geo2gsm", "GEO_to_GSM", "geo2j2000", ...; names may contain '2', so match known names
function parse_coord_transform(s::AbstractString)
    u = uppercase(replace(s, r"_to_"i => "2"))
    for c in COORD_SYSTEMS
        prefix = string(c, "2")
        startswith(u, prefix) || continue
        out = Symbol(chopprefix(u, prefix))
        out in COORD_SYSTEMS && return c, out
    end
    throw(ArgumentError("Could not parse coordinate system conversion string: '$s'. Expected format like 'geo2gsm'."))
end

function clean_posit!(posit::Array{T, 3}, Blocal::Matrix{T}, Nposit) where {T}
    for (i, n) in enumerate(Nposit)
        posit[:, (n + 1):end, i] .= T(NaN)
        Blocal[(n + 1):end, i] .= T(NaN)
    end
    return
end
