"""
    get_datetime(X::Dict)

Extract datetime from input dictionary X.
Supports 'dateTime' or 'Time' keys with DateTime or String values.
"""
function get_datetime(X::Dict)
    if haskey(X, "dateTime")
        dt_val = X["dateTime"]
    elseif haskey(X, "Time")
        dt_val = X["Time"]
    else
        error("No date/time information found in input dictionary. Expected 'dateTime' or 'Time' key.")
    end

    if isa(dt_val, DateTime)
        return dt_val
    elseif isa(dt_val, String)
        return DateTime(dt_val)
    else
        error("Date/time value must be a DateTime or String, got $(typeof(dt_val))")
    end
end

"""
    process_coords_time(X::Dict)

Process coordinates and time from input dictionary X.
Returns ntime, iyear, idoy, ut, x1, x2, x3 arrays for IRBEM functions.
"""
function process_coords_time(X::Dict)
    # Handle array or scalar inputs
    if any(isa(X[k], AbstractArray) for k in ["x1", "x2", "x3"] if haskey(X, k))
        # Array input
        if haskey(X, "dateTime") && isa(X["dateTime"], AbstractArray)
            ntime = length(X["dateTime"])
        elseif haskey(X, "Time") && isa(X["Time"], AbstractArray)
            ntime = length(X["Time"])
        else
            ntime = length(X["x1"])
        end

        # Initialize arrays
        iyear = zeros(Int32, ntime)
        idoy = zeros(Int32, ntime)
        ut = zeros(Float64, ntime)
        x1 = zeros(Float64, ntime)
        x2 = zeros(Float64, ntime)
        x3 = zeros(Float64, ntime)

        # Process each time point
        for i in 1:ntime
            # Extract datetime
            if haskey(X, "dateTime")
                dt_val = isa(X["dateTime"], AbstractArray) ? X["dateTime"][i] : X["dateTime"]
            else
                dt_val = isa(X["Time"], AbstractArray) ? X["Time"][i] : X["Time"]
            end

            # Convert to DateTime if needed
            dt = isa(dt_val, DateTime) ? dt_val : DateTime(dt_val)

            # Set time values
            iyear[i] = year(dt)
            idoy[i] = dayofyear(dt)
            ut[i] = hour(dt) * 3600 + minute(dt) * 60 + second(dt)

            # Set coordinate values
            x1[i] = isa(X["x1"], AbstractArray) ? X["x1"][i] : X["x1"]
            x2[i] = isa(X["x2"], AbstractArray) ? X["x2"][i] : X["x2"]
            x3[i] = isa(X["x3"], AbstractArray) ? X["x3"][i] : X["x3"]
        end
    else
        # Scalar input
        ntime = 1

        # Get datetime
        dt = get_datetime(X)

        # Initialize arrays
        iyear = [Int32(year(dt))]
        idoy = [Int32(dayofyear(dt))]
        ut = [Float64(hour(dt) * 3600 + minute(dt) * 60 + second(dt))]
        x1 = [Float64(X["x1"])]
        x2 = [Float64(X["x2"])]
        x3 = [Float64(X["x3"])]
    end

    return Int32(ntime), iyear, idoy, ut, x1, x2, x3
end

"""
    prepare_maginput(maginput::Dict, ntime::Int)

Process magnetic field model inputs from input dictionary.
Returns a properly formatted array for IRBEM functions.
"""
function prepare_maginput(maginput::Dict, ntime)
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
    lookup_table = Dict(
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
function clean_posit!(posit, Nposit)
    for i in eachindex(Nposit)
        n = Nposit[i]
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