_vec(x) = [x]
_vec(x::AbstractVector) = x
_only(x) = length(x) == 1 ? x[1] : x
_only(A::AbstractMatrix) = size(A, 2) == 1 ? A[:, 1] : A

"""
    get_datetime(X::Dict)

Extract datetime from input dictionary X.
Supports 'dateTime' or 'Time' keys with DateTime or String values.
"""
function get_datetime(X::Dict)
    # Extract date/time value (scalar or array)
    dt_val = haskey(X, "dateTime") ? X["dateTime"] : haskey(X, "Time") ? X["Time"] : nothing
    if dt_val === nothing
        error("No date/time information found in input dictionary. Expected 'dateTime' or 'Time' key.")
    end
    # Handle array or scalar
    if isa(dt_val, AbstractArray)
        return [isa(dt, DateTime) ? dt : DateTime(dt) for dt in dt_val]
    elseif isa(dt_val, DateTime)
        return dt_val
    elseif isa(dt_val, String)
        return DateTime(dt_val)
    else
        error("Date/time value must be a DateTime, String, or Array thereof, got $(typeof(dt_val))")
    end
end

prepare_irbem(time, x, coord="GDZ", maginput=Dict(); kext=KEXT[], options=OPTIONS[]) = (
    ntime(time), parse_kext(kext), options, coord_sys(coord),
    decompose_time(time)..., prepare_loc(x)...,
    prepare_maginput(maginput)
)

function prepare_irbem(model::MagneticField, X, maginput=Dict())
    time = get_datetime(X)
    return (
        ntime(time), model.kext, model.options, model.sysaxes,
        decompose_time(time)..., prepare_loc(X)...,
        prepare_maginput(maginput)
    )
end

function decompose_time(x::AbstractVector)
    dt = DateTime.(x)
    iyear = Int32.(year.(dt))
    idoy = Int32.(dayofyear.(dt))
    ut = Float64.(hour.(dt) * 3600 + minute.(dt) * 60 + second.(dt) + millisecond.(dt) / 1000)
    iyear, idoy, ut
end

decompose_time(dt) = decompose_time(_vec(dt))

function prepare_time(dt::AbstractVector)
    ntime = Int32(length(dt))
    ntime, decompose_time(dt)...
end

prepare_time(dt::DateTime) = prepare_time([dt])

ntime(time) = Int32(1)
ntime(time::AbstractVector) = Int32(length(time))

function prepare_loc(x1, x2, x3)
    _x1 = x1 isa Number ? [Float64(x1)] : convert(Vector{Float64}, x1)
    _x2 = x2 isa Number ? [Float64(x2)] : convert(Vector{Float64}, x2)
    _x3 = x3 isa Number ? [Float64(x3)] : convert(Vector{Float64}, x3)
    _x1, _x2, _x3
end

prepare_loc(x::AbstractArray) = prepare_loc(x[1, :], x[2, :], x[3, :])
prepare_loc(x::AbstractVector{<:AbstractVector}) = prepare_loc([getindex.(x, i) for i in 1:3]...)
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
    prepare_maginput(maginput::Dict, ntime::Int)

Process magnetic field model inputs from input dictionary.
Returns a properly formatted array for IRBEM functions.
"""
function prepare_maginput(maginput::Dict, ntime=nothing)
    # IRBEM expects a 25-element array for maginput
    maginput_array = zeros(Float64, 25)

    # Map of parameter names to indices in the maginput array
    param_indices = Dict(
        "Kp" => 1,
        "Dst" => 2,
        "dens" => 3,
        "velo" => 4,
        "Pdyn" => 5,
        "ByIMF" => 6,
        "BzIMF" => 7,
        "G1" => 8,
        "G2" => 9,
        "G3" => 10,
        "W1" => 11,
        "W2" => 12,
        "W3" => 13,
        "W4" => 14,
        "W5" => 15,
        "W6" => 16,
        "AL" => 17
    )

    # Fill the array with values from the input dictionary
    for (param, idx) in param_indices
        if haskey(maginput, param)
            val = maginput[param]
            if isa(val, AbstractArray)
                # If array input, use the first value
                maginput_array[idx] = Float64(val[1])
            else
                # Otherwise use the scalar value
                maginput_array[idx] = Float64(val)
            end
        end
    end

    return maginput_array
end


parse_kext(kext::Integer) = kext
function parse_kext(kext)
    idx = findfirst(isequal(kext), EXT_MODELS)
    !isnothing(idx) ? idx - 1 : throw(ArgumentError("Unknown external field model: $kext. Valid models are $EXT_MODELS"))
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
    lookup_table = Dict{String,Int32}(
        "GDZ" => 0,
        "GEO" => 1,
        "GSM" => 2,
        "GSE" => 3,
        "SM" => 4,
        "GEI" => 5,
        "MAG" => 6,
        "SPH" => 7,
        "RLL" => 8
    )
    get(lookup_table, uppercase(axes)) do
        error("Unknown coordinate system: $axes. Choose from GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL.")
    end
end

coord_sys(x::Integer) = Int32(x)

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
beta(Ek, Erest=511.0) = sqrt(1 - ((Ek / Erest) + 1)^(-2))

"""
    gamma(Ek, Erest=511.0)

Calculate relativistic gamma factor for a particle with kinetic energy Ek.
Ek and Erest must be in the same units (default is keV).
"""
gamma(Ek, Erest=511.0) = 1 / sqrt(1 - beta(Ek, Erest)^2)

"""
    vparallel(Ek, Bm, B, Erest=511.0)

Calculate parallel velocity for a particle with kinetic energy Ek,
at a location with magnetic field B, with mirror point field Bm.
Ek and Erest must be in the same units (default is keV).
Returns velocity in m/s.
"""
vparallel(Ek, Bm, B, Erest=511.0) = c * beta(Ek, Erest) * sqrt(1 - abs(B / Bm))

"""
    clean_posit!(posit, Nposit)

Remove trailing NaN values from the posit array.
"""
function clean_posit!(posit, Nposit::AbstractVector)
    for (i, n) in enumerate(Nposit)
        posit[:, n+1:end, i] .= NaN
    end
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
    quote
        $(Expr(:block, [:($(esc(v)) = Ref{$T}()) for v in vars]...))
    end
end