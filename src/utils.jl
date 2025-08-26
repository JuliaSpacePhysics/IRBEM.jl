const SF64 = Scalar{Float64}
const RF64 = Ref{Float64}

# https://docs.julialang.org/en/v1/base/c/#Base.unsafe_convert
# We should make sure the input is a c-contiguous Float64 array
@inline vecf(x::Vector{Float64}) = x
@inline vecf(x) = convert(Vector{Float64}, x)
@inline vecf(x::StaticArray) = eltype(x) == Float64 ? x : Float64.(x)
vecf(x::Number) = SF64(x)
@inline arrf(x) = eltype(x) == Float64 ? x : convert(Array{Float64}, x)
_vec(x) = [x]
_vec(x::AbstractVector) = x
_first(x::AbstractVector) = x[1]
_first(A::AbstractMatrix) = A[:, 1]
_first(A::AbstractArray{T, 3}) where {T} = A[:, :, 1]

_deref(x::Ref) = x[]
_deref(x) = x

"""
    get_datetime(X::Dict)

Extract datetime from input dictionary X.
Supports 'dateTime' or 'Time' keys with DateTime or String values.
"""
function get_datetime(X::Dict)
    dt_val = get(X, "dateTime", get(X, "Time", nothing))
    return !isnothing(dt_val) ? Dates.DateTime.(dt_val) : error("No date/time information found in input dictionary. Expected 'dateTime' or 'Time' key.")
end

prepare_irbem(time, x, coord = "GDZ", maginput = Dict(); kext = KEXT[], options = OPTIONS[]) = (
    ntime(time), parse_kext(kext), options, coord_sys(coord),
    decompose_time(time)..., prepare_loc(x)...,
    prepare_maginput(maginput),
)

prepare_irbem(time, x::CoordinateVector, maginput = (;); kext = KEXT[], options = OPTIONS[]) = (
    ntime(time), parse_kext(kext), options, coord_sys(x),
    decompose_time(time)..., prepare_loc(x)...,
    prepare_maginput(maginput),
)

function prepare_irbem(model::MagneticField, X::AbstractDict, maginput = Dict())
    time = get_datetime(X)
    return (
        ntime(time), model.kext, model.options, model.sysaxes,
        decompose_time(time)..., prepare_loc(X)...,
        prepare_maginput(maginput),
    )
end

"""
    decompose_time_s(dt::DateTime)

Decompose a single DateTime into year, day of year, and UT.
"""
function decompose_time_s(dt::DateTime)
    iyear = Int32(year(dt))
    idoy = Int32(dayofyear(dt))
    ut = Float64(hour(dt) * 3600 + minute(dt) * 60 + second(dt) + millisecond(dt) / 1000)
    return iyear, idoy, ut
end

decompose_time_s(dt) = decompose_time_s(DateTime(dt))

function decompose_time(x::AbstractVector)
    dt = eltype(x) <: DateTime ? x : DateTime.(x)
    iyear = @. Int32(year(dt))
    idoy = @. Int32(dayofyear(dt))
    ut = @. Float64(hour(dt) * 3600 + minute(dt) * 60 + second(dt) + millisecond(dt) / 1000)
    return iyear, idoy, ut
end

decompose_time(dt) = decompose_time(_vec(dt))

function prepare_time(dt::AbstractVector)
    ntime = Int32(length(dt))
    return ntime, decompose_time(dt)...
end

prepare_time(dt::DateTime) = prepare_time([dt])

ntime(time) = Int32(1)
ntime(time::AbstractVector) = Int32(length(time))

prepare_loc(x1, x2, x3) = vecf(x1), vecf(x2), vecf(x3)
prepare_loc(x::AbstractVector) = SF64(x[1]), SF64(x[2]), SF64(x[3])
prepare_loc(x::AbstractArray) = vecf(x[1, :]), vecf(x[2, :]), vecf(x[3, :])
prepare_loc(x::AbstractVector{<:AbstractVector}) = prepare_loc((getindex.(x, i) for i in 1:3)...)
prepare_loc(X::Dict) = prepare_loc(X["x1"], X["x2"], X["x3"])


"""
    process_coords_time(X::Dict)

Process coordinates and time from input dictionary X.
Returns ntime, iyear, idoy, ut, x1, x2, x3 arrays for IRBEM functions.
"""
function process_coords_time(X::Dict)
    ntime, iyear, idoy, ut = prepare_time(get_datetime(X))
    x1, x2, x3 = prepare_loc(X)
    return ntime, iyear, idoy, ut, x1, x2, x3
end

"""
    with_case_variants(dict::Dict{S, V}) where {S <: AbstractString, V}

Creates a new dictionary with both uppercase and lowercase variants of each key.
The original keys and values are preserved, and lowercase/uppercase variants are added.

Example:
```julia
original = Dict("ABC" => 1, "DEF" => 2)
result = with_case_variants(original)
# result has keys: "ABC", "abc", "DEF", "def"
```
"""
function with_case_variants(dict)
    result = copy(dict)
    for (key, value) in dict
        result[lowercase(key)] = value
        result[uppercase(key)] = value
    end
    return result
end

_get_param(maginput::Dict{K, V}, param, default = nothing) where {K, V} = get(maginput, K(param), default)
_get_param(maginput, param, default = nothing) = get(maginput, param, default)

"""
    prepare_maginput(maginput)

Process magnetic field model inputs from input dictionary.
Returns a properly formatted array for IRBEM functions.
"""
function prepare_maginput(maginput)
    # IRBEM expects a 25-element array for maginput
    out = MVector{25, Float64}(undef)
    # Fill the array with values from the input dictionary
    for (idx, param) in enumerate(param_indices)
        val = _get_param(maginput, param, 0)
        # If array input, use the first value
        out[idx] = Float64(first(val))
    end

    return out
end


parse_kext(kext::Integer) = kext
function parse_kext(kext)
    idx = findfirst(isequal(kext), EXT_MODELS)
    return !isnothing(idx) ? idx - 1 : throw(ArgumentError("Unknown external field model: $kext. Valid models are $EXT_MODELS"))
end

"""
    coord_sys(axes)

Look up the IRBEM coordinate system integer given a string representation.

Coordinate systems:
- 0: GDZ: (altitude, latitude, east longitude - km, deg, deg)
- 1: GEO: Cartesian GEO - Re
- 2: GSM: Cartesian GSM - Re
- 3: GSE: Cartesian GSE - Re
- 4: SM: Cartesian SM - Re
- 5: GEI: Cartesian GEI - Re
- 6: MAG: Cartesian MAG - Re
- 7: SPH: Spherical GEO - (radial distance, latitude, east longitude - Re, deg, deg)
- 8: RLL: Spherical GEO - (radial distance, latitude, east longitude - Re, deg, deg)

Returns the corresponding integer code.
"""
function coord_sys(axes)
    return get(coord_sys_lookup, axes) do
        error("Unknown coordinate system: $axes. Choose from GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL.")
    end
end

function coord_sys(x::Symbol)
    return get(coord_sys_lookup_sym, x) do
        error("Unknown coordinate system: $x. Choose from GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL.")
    end
end

coord_sys(x::Integer) = Int32(x)
coord_sys(::Type{S}) where {S <: AbstractCoordinateSystem} = coord_sys(nameof(S))
coord_sys(x::AbstractCoordinateSystem) = coord_sys(typeof(x))

parse_coord_transform(pair) = pair[1], pair[2]
function parse_coord_transform(s::String)
    # Accept formats like "geo2gsm", "GEO2GSM", "geo_to_gsm", etc.
    s_clean = replace(s, "_to_" => "2")
    if occursin("2", s_clean)
        parts = split(s_clean, "2")
        if length(parts) == 2
            return parts[1], parts[2]
        end
    end
    error("Could not parse coordinate system conversion string: '$s'. Expected format like 'geo2gsm'.")
end

# Helper functions for relativistic calculations
"""
    beta(Ek, Erest=511.0)

Calculate relativistic beta (v/c) for a particle with kinetic energy Ek.
Ek and Erest must be in the same units (default is keV).
"""
beta(Ek, Erest = 511.0) = sqrt(1 - ((Ek / Erest) + 1)^(-2))

"""
    gamma(Ek, Erest=511.0)

Calculate relativistic gamma factor for a particle with kinetic energy Ek.
Ek and Erest must be in the same units (default is keV).
"""
gamma(Ek, Erest = 511.0) = 1 / sqrt(1 - beta(Ek, Erest)^2)

"""
    vparallel(Ek, Bm, B, Erest=511.0)

Calculate parallel velocity for a particle with kinetic energy Ek,
at a location with magnetic field B, with mirror point field Bm.
Ek and Erest must be in the same units (default is keV).
Returns velocity in m/s.
"""
vparallel(Ek, Bm, B, Erest = 511.0) = c * beta(Ek, Erest) * sqrt(1 - abs(B / Bm))

"""
    clean_posit!(posit, Nposit)

Remove trailing NaN values from the posit array.
"""
function clean_posit!(posit, Nposit::AbstractVector)
    for (i, n) in enumerate(Nposit)
        posit[:, (n + 1):end, i] .= NaN
    end
    return
end

"""
    @init_refs(Type, var1, var2, ...)

Creates variables var1, var2, ... each initialized as `Ref{Type}()`.
Example:
    @init_refs(Float64, 0.0, Lm, Lstar)
expands to:
    Lm = Ref{Float64}(0.0)
    Lstar = Ref{Float64}(0.0)
"""
macro init_refs(T, vars...)
    return quote
        $(Expr(:block, [:($(esc(v)) = Ref{$T}()) for v in vars]...))
    end
end
